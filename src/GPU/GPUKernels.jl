"""
GPU kernel implementations for tensor-product element stiffness matvec.

Uses KernelAbstractions.jl for GPU-portable kernel definitions.
Each kernel processes one element per workgroup; threads within the workgroup
correspond to the (i,j) nodes of the 2D element (nnodes×nnodes layout).

Design:
- Group index = element
- Local index = node within element (linearized ij, 1..nnodes²)
- Metric weights g[:,:,:,el] are read from global memory (element-contiguous ✓)
- H, Ht are small matrices (5×5 for p=4) — loaded once per workgroup
- Scatter-add uses @atomic since boundary DOFs are shared between elements

Memory layout for g:
  Antiplane:   g[i,j,1..3,el]
  Plane-strain: g[i,j,1..10,el]

The tensor-product formula for antiplane (and each block of plane-strain) is:
  f = H*(Wξξ.*(Hᵀ*u)) + (Wηη.*(u*H))*Hᵀ
    + H*(Wξη.*(u*H)) + (Wξη.*(Hᵀ*u))*Hᵀ

This kernel evaluates the same formula but distributes the nnodes² nodes across
GPU threads, using shared memory for the element-local u_el and intermediate results.
"""

using KernelAbstractions
using KernelAbstractions.Extras: @unroll


# ============================================================================
# Antiplane kernel
# ============================================================================

"""
    apply_stiffness_antiplane_kernel!(y, x, g, H, Ht, dof_id, nnodes)

KernelAbstractions kernel: one workgroup per element, nnodes² threads per workgroup.

Computes y += K_el * x_el for each element via tensor-product matvec, then
scatter-adds (atomic) the result into the global y vector.

# Arguments
- `y`: Global output vector (length ndof)
- `x`: Global input vector (length ndof)
- `g`: Metric weights [nnodes, nnodes, 3, n_elements]
- `H`: Derivative matrix [nnodes, nnodes]
- `Ht`: H' [nnodes, nnodes]
- `dof_id`: Connectivity [nnodes, nnodes, n_elements]
- `nnodes`: Nodes per direction (polynomial_degree + 1)
"""
@kernel function apply_stiffness_antiplane_kernel!(y, x, g, H, Ht, dof_id, nnodes)
    el = @index(Group, Linear)
    ij = @index(Local, Linear)

    # Compute 2D (i,j) index from linear thread index
    j = cld(ij, nnodes)         # column (η direction)
    i = ij - (j - 1) * nnodes   # row (ξ direction)

    # Shared memory for element-local arrays
    u_el  = @localmem eltype(x) (nnodes * nnodes,)
    f_el  = @localmem eltype(y) (nnodes * nnodes,)
    tmp1  = @localmem eltype(x) (nnodes * nnodes,)
    tmp2  = @localmem eltype(x) (nnodes * nnodes,)

    # Load u_el from global x via connectivity
    @inbounds u_el[ij] = x[dof_id[i, j, el]]
    f_el[ij] = zero(eltype(y))
    @synchronize

    T = eltype(x)

    # Tensor-product matvec: 4 terms
    # Term 1: H * (Wξξ .* (Hᵀ * u_el))
    # Term 2: (Wηη .* (u_el * H)) * Hᵀ
    # Term 3: H * (Wξη .* (u_el * H))    [cross term 1]
    # Term 4: (Wξη .* (Hᵀ * u_el)) * Hᵀ [cross term 2]

    # --- dU/dξ = Hᵀ * u_el ---
    # (Hᵀ * u_el)[i,j] = Σ_k Ht[i,k] * u_el[k,j] = Σ_k H[k,i] * u_el[k,j]
    val_dxi = zero(T)
    @unroll for k in 1:nnodes
        val_dxi += Ht[i, k] * u_el[(j-1)*nnodes + k]
    end
    tmp1[ij] = val_dxi
    @synchronize

    # --- dU/dη = u_el * H ---
    # (u_el * H)[i,j] = Σ_k u_el[i,k] * H[k,j]
    val_deta = zero(T)
    @unroll for k in 1:nnodes
        val_deta += u_el[(k-1)*nnodes + i] * H[k, j]
    end
    tmp2[ij] = val_deta
    @synchronize

    # Scale by metric weights
    @inbounds begin
        Wxi  = g[i, j, 1, el]   # Wξξ
        Weta = g[i, j, 2, el]   # Wηη
        Wxie = g[i, j, 3, el]   # Wξη = Wηξ

        # tmp1 *= Wξξ,  tmp2 *= Wηη  (reuse arrays for first pair)
        tmp1[ij] = tmp1[ij] * Wxi
        tmp2[ij] = tmp2[ij] * Weta
    end
    @synchronize

    # f += H * tmp1  (H * Wξξ*dU/dξ)
    val_f = zero(T)
    @unroll for k in 1:nnodes
        val_f += H[i, k] * tmp1[(j-1)*nnodes + k]
    end
    # f += tmp2 * Hᵀ  (Wηη*dU/dη * Hᵀ)
    @unroll for k in 1:nnodes
        val_f += tmp2[(k-1)*nnodes + i] * Ht[k, j]
    end
    f_el[ij] = val_f
    @synchronize

    # Cross-term 1: H * (Wξη * dU/dη)
    # Restore dU/dη into tmp2 (already done above in val_deta, but tmp2 was overwritten)
    val_deta2 = zero(T)
    @unroll for k in 1:nnodes
        val_deta2 += u_el[(k-1)*nnodes + i] * H[k, j]
    end
    tmp2[ij] = val_deta2
    @synchronize

    @inbounds tmp2[ij] = tmp2[ij] * g[i, j, 3, el]  # Wξη * dU/dη
    @synchronize

    val_cross1 = zero(T)
    @unroll for k in 1:nnodes
        val_cross1 += H[i, k] * tmp2[(j-1)*nnodes + k]
    end

    # Cross-term 2: (Wξη * dU/dξ) * Hᵀ
    # Restore dU/dξ into tmp1
    val_dxi2 = zero(T)
    @unroll for k in 1:nnodes
        val_dxi2 += Ht[i, k] * u_el[(j-1)*nnodes + k]
    end
    tmp1[ij] = val_dxi2
    @synchronize

    @inbounds tmp1[ij] = tmp1[ij] * g[i, j, 3, el]  # Wξη * dU/dξ
    @synchronize

    val_cross2 = zero(T)
    @unroll for k in 1:nnodes
        val_cross2 += tmp1[(k-1)*nnodes + i] * Ht[k, j]
    end

    f_el[ij] += val_cross1 + val_cross2
    @synchronize

    # Scatter-add into global y (atomic: multiple elements share boundary DOFs)
    @inbounds CUDA.@atomic y[dof_id[i, j, el]] += f_el[ij]
end


# ============================================================================
# Antiplane CPU-side launch function
# ============================================================================

"""
    launch_stiffness_antiplane!(y, x, g, H, Ht, dof_id, n_elements, nnodes)

Launch the antiplane tensor-product stiffness kernel on GPU.

Zeroes y first, then launches one workgroup per element with nnodes² threads.
"""
function launch_stiffness_antiplane!(
    y::AbstractVector{T},
    x::AbstractVector{T},
    g::AbstractArray{T,4},
    H::AbstractMatrix{T},
    Ht::AbstractMatrix{T},
    dof_id::AbstractArray{<:Integer,3},
    n_elements::Int,
    nnodes::Int
) where T
    fill!(y, zero(T))
    kernel = apply_stiffness_antiplane_kernel!(CUDABackend())
    kernel(y, x, g, H, Ht, dof_id, nnodes;
           ndrange = (n_elements * nnodes * nnodes,),
           groupsize = nnodes * nnodes)
    KernelAbstractions.synchronize(CUDABackend())
    return y
end


# ============================================================================
# Plane-strain kernel (K_xx, K_yy, K_xy blocks)
# ============================================================================

"""
    apply_stiffness_ps_block_kernel!(f, u, g, H, Ht, dof_id, nnodes, c1, c2, c3, c4, mode)

KernelAbstractions kernel for one block contribution in the plane-strain matvec.

mode=0: symmetric block (K_xx or K_yy): uses components c1(ξξ), c2(ηη), c3(ξη), c4=c3(ηξ)
mode=1: asymmetric block (K_xy): c1..c4 are all distinct
mode=2: transposed asymmetric block (K_yx = K_xy^T): c3 and c4 swapped vs K_xy
"""
@kernel function apply_stiffness_ps_block_kernel!(
    f, u, g, H, Ht, dof_id, nnodes, c1, c2, c3, c4, ndof_offset_u, ndof_offset_f
)
    el = @index(Group, Linear)
    ij = @index(Local, Linear)

    j = cld(ij, nnodes)
    i = ij - (j - 1) * nnodes

    u_el  = @localmem eltype(u) (nnodes * nnodes,)
    f_el  = @localmem eltype(f) (nnodes * nnodes,)
    tmp1  = @localmem eltype(u) (nnodes * nnodes,)
    tmp2  = @localmem eltype(u) (nnodes * nnodes,)

    T = eltype(u)

    # Load u_el using spatial dof_id (offset by ndof for y-component)
    @inbounds begin
        idx = dof_id[i, j, el]
        u_el[ij] = u[ndof_offset_u + idx]
    end
    f_el[ij] = zero(T)
    @synchronize

    # dU/dξ = Hᵀ * u_el
    val_dxi = zero(T)
    @unroll for k in 1:nnodes
        val_dxi += Ht[i, k] * u_el[(j-1)*nnodes + k]
    end
    tmp1[ij] = val_dxi

    # dU/dη = u_el * H
    val_deta = zero(T)
    @unroll for k in 1:nnodes
        val_deta += u_el[(k-1)*nnodes + i] * H[k, j]
    end
    tmp2[ij] = val_deta
    @synchronize

    @inbounds begin
        tmp1[ij] = tmp1[ij] * g[i, j, c1, el]   # c1 * dU/dξ
        tmp2[ij] = tmp2[ij] * g[i, j, c2, el]   # c2 * dU/dη
    end
    @synchronize

    # H * (c1*dU/dξ)
    val_f = zero(T)
    @unroll for k in 1:nnodes
        val_f += H[i, k] * tmp1[(j-1)*nnodes + k]
    end
    # (c2*dU/dη) * Hᵀ
    @unroll for k in 1:nnodes
        val_f += tmp2[(k-1)*nnodes + i] * Ht[k, j]
    end
    f_el[ij] = val_f
    @synchronize

    # Cross-term 1: H * (c3 * dU/dη)
    val_deta2 = zero(T)
    @unroll for k in 1:nnodes
        val_deta2 += u_el[(k-1)*nnodes + i] * H[k, j]
    end
    tmp2[ij] = val_deta2
    @synchronize

    @inbounds tmp2[ij] = tmp2[ij] * g[i, j, c3, el]
    @synchronize

    val_cross1 = zero(T)
    @unroll for k in 1:nnodes
        val_cross1 += H[i, k] * tmp2[(j-1)*nnodes + k]
    end

    # Cross-term 2: (c4 * dU/dξ) * Hᵀ
    val_dxi2 = zero(T)
    @unroll for k in 1:nnodes
        val_dxi2 += Ht[i, k] * u_el[(j-1)*nnodes + k]
    end
    tmp1[ij] = val_dxi2
    @synchronize

    @inbounds tmp1[ij] = tmp1[ij] * g[i, j, c4, el]
    @synchronize

    val_cross2 = zero(T)
    @unroll for k in 1:nnodes
        val_cross2 += tmp1[(k-1)*nnodes + i] * Ht[k, j]
    end

    f_el[ij] += val_cross1 + val_cross2
    @synchronize

    @inbounds CUDA.@atomic f[ndof_offset_f + dof_id[i, j, el]] += f_el[ij]
end


"""
    launch_stiffness_plane_strain!(y, x, g, H, Ht, dof_id, n_elements, nnodes, ndof)

Launch the plane-strain tensor-product stiffness kernel on GPU.

Zeroes y first, then launches 4 kernel passes for K_xx, K_yy, K_xy, K_yx.
"""
function launch_stiffness_plane_strain!(
    y::AbstractVector{T},
    x::AbstractVector{T},
    g::AbstractArray{T,4},
    H::AbstractMatrix{T},
    Ht::AbstractMatrix{T},
    dof_id::AbstractArray{<:Integer,3},
    n_elements::Int,
    nnodes::Int,
    ndof::Int
) where T
    fill!(y, zero(T))
    groupsize = nnodes * nnodes
    ndrange = (n_elements * groupsize,)
    backend = CUDABackend()

    kernel = apply_stiffness_ps_block_kernel!(backend)

    # K_xx * ux → fx  (components 1,2,3,3; offset_u=0, offset_f=0)
    kernel(y, x, g, H, Ht, dof_id, nnodes, 1, 2, 3, 3, 0, 0;
           ndrange=ndrange, groupsize=groupsize)
    KernelAbstractions.synchronize(backend)

    # K_yy * uy → fy  (components 4,5,6,6; offset_u=ndof, offset_f=ndof)
    kernel(y, x, g, H, Ht, dof_id, nnodes, 4, 5, 6, 6, ndof, ndof;
           ndrange=ndrange, groupsize=groupsize)
    KernelAbstractions.synchronize(backend)

    # K_xy * uy → fx  (components 7,8,9,10; offset_u=ndof, offset_f=0)
    kernel(y, x, g, H, Ht, dof_id, nnodes, 7, 8, 9, 10, ndof, 0;
           ndrange=ndrange, groupsize=groupsize)
    KernelAbstractions.synchronize(backend)

    # K_yx * ux → fy  (= K_xy^T: swap c9↔c10; offset_u=0, offset_f=ndof)
    kernel(y, x, g, H, Ht, dof_id, nnodes, 7, 8, 10, 9, 0, ndof;
           ndrange=ndrange, groupsize=groupsize)
    KernelAbstractions.synchronize(backend)

    return y
end
