"""
Matrix-free stiffness operator for iterative solvers.

Implements efficient element-by-element matrix-vector products without
assembling the global stiffness matrix. Uses tensor-product formulation
with compact metric weight arrays (MetricWeightsAntiplane/PlaneStrain)
instead of stored full element stiffness matrices.

Memory reduction vs old K_el approach:
- Antiplane p=4: ~40× smaller (75 doubles vs 3125 doubles per element)
- Plane-strain p=4: ~160× smaller
"""

"""
    apply_stiffness!(y, x, weights, H, Ht, dof_id, n_elements)

Apply stiffness operator y = K*x using element-by-element tensor products (matrix-free).

# Arguments
- `y::Vector`: Output vector (modified in-place, zeroed first)
- `x::Vector`: Input vector
- `weights::MetricWeightsAntiplane{T}`: Compact metric weight arrays
- `H::Matrix{T}`: Derivative matrix [nnodes, nnodes]
- `Ht::Matrix{T}`: H' (pre-transposed)
- `dof_id::Array{Int,3}`: Connectivity matrix [nnodes, nnodes, n_elements]
- `n_elements::Int`: Number of elements

# Algorithm
```
y .= 0
for each element:
    gather u_el from global x
    compute f_el via tensor-product: apply_element_stiffness_antiplane!
    scatter f_el into global y
```
"""
function apply_stiffness!(
    y::Vector{T},
    x::Vector{T},
    weights::MetricWeightsAntiplane{T},
    H::Matrix{T},
    Ht::Matrix{T},
    dof_id::Array{Int,3},
    n_elements::Int
) where {T<:AbstractFloat}
    fill!(y, zero(T))

    nnodes = weights.nnodes
    u_el  = zeros(T, nnodes, nnodes)
    f_el  = zeros(T, nnodes, nnodes)
    tmp1  = zeros(T, nnodes, nnodes)
    tmp2  = zeros(T, nnodes, nnodes)

    for el in 1:n_elements
        local_idx = dof_id[:, :, el]

        # Gather: fill u_el from global x
        @inbounds for j in 1:nnodes, i in 1:nnodes
            u_el[i, j] = x[local_idx[i, j]]
        end

        # Tensor-product element matvec
        fill!(f_el, zero(T))
        apply_element_stiffness_antiplane!(f_el, u_el, view(weights.g, :, :, :, el),
                                           H, Ht, tmp1, tmp2)

        # Scatter-add into global y
        @inbounds for j in 1:nnodes, i in 1:nnodes
            y[local_idx[i, j]] += f_el[i, j]
        end
    end

    return y
end


"""
    StiffnessOperator{T} <: AbstractMatrix{T}

Matrix-free stiffness operator for use with iterative linear solvers.

Uses compact metric weight arrays and tensor-product matvec instead of stored
full element stiffness matrices. This enables GPU-ready element loops.

# Fields
- `weights::MetricWeightsAntiplane{T}`: Compact metric weights (replaces K_el)
- `H::Matrix{T}`: Derivative matrix [nnodes, nnodes]
- `Ht::Matrix{T}`: H' (pre-transposed for efficient mul!)
- `dof_id::Array{Int,3}`: Connectivity matrix [nnodes, nnodes, n_elements]
- `n_elements::Int`: Number of elements
- `ndof::Int`: Total number of global DOFs
- `fltni::Vector{Int}`: Non-fault/creep DOF indices (interior DOFs)
- `x_full::Vector{T}`: Workspace for full DOF vector (pre-allocated)
- `y_full::Vector{T}`: Workspace for full DOF vector (pre-allocated)
- `u_el, f_el, tmp1, tmp2::Matrix{T}`: Element-local workspaces [nnodes, nnodes]
"""
struct StiffnessOperator{T<:AbstractFloat} <: AbstractMatrix{T}
    weights::MetricWeightsAntiplane{T}
    H::Matrix{T}
    Ht::Matrix{T}
    dof_id::Array{Int,3}
    n_elements::Int
    ndof::Int
    fltni::Vector{Int}
    x_full::Vector{T}
    y_full::Vector{T}
    # Element-local workspaces (avoid per-element allocation in mul!)
    u_el::Matrix{T}
    f_el::Matrix{T}
    tmp1::Matrix{T}
    tmp2::Matrix{T}
end

"""
    StiffnessOperator(weights, H, Ht, dof_id, n_elements, ndof, fltni)

Construct matrix-free stiffness operator with pre-allocated workspaces.
"""
function StiffnessOperator(
    weights::MetricWeightsAntiplane{T},
    H::Matrix{T},
    Ht::Matrix{T},
    dof_id::Array{Int,3},
    n_elements::Int,
    ndof::Int,
    fltni::Vector{Int}
) where {T<:AbstractFloat}
    x_full = zeros(T, ndof)
    y_full = zeros(T, ndof)
    nnodes = weights.nnodes
    u_el = zeros(T, nnodes, nnodes)
    f_el = zeros(T, nnodes, nnodes)
    tmp1 = zeros(T, nnodes, nnodes)
    tmp2 = zeros(T, nnodes, nnodes)
    return StiffnessOperator(weights, H, Ht, dof_id, n_elements, ndof, fltni,
                             x_full, y_full, u_el, f_el, tmp1, tmp2)
end

"""
    LinearAlgebra.mul!(y, A::StiffnessOperator, x)

Matrix-vector product y = A*x for matrix-free stiffness operator (antiplane).

1. Expand x to full DOF vector (zero on fault/creep boundaries)
2. Apply stiffness element-by-element via tensor products
3. Extract non-fault DOFs from result
"""
function LinearAlgebra.mul!(y::Vector, A::StiffnessOperator{T}, x::Vector) where {T}
    fill!(A.x_full, zero(T))
    A.x_full[A.fltni] .= x

    fill!(A.y_full, zero(T))

    nnodes = A.weights.nnodes

    for el in 1:A.n_elements
        local_idx = A.dof_id[:, :, el]

        # Gather
        @inbounds for j in 1:nnodes, i in 1:nnodes
            A.u_el[i, j] = A.x_full[local_idx[i, j]]
        end

        # Tensor-product matvec
        fill!(A.f_el, zero(T))
        apply_element_stiffness_antiplane!(A.f_el, A.u_el,
                                           view(A.weights.g, :, :, :, el),
                                           A.H, A.Ht, A.tmp1, A.tmp2)

        # Scatter-add
        @inbounds for j in 1:nnodes, i in 1:nnodes
            A.y_full[local_idx[i, j]] += A.f_el[i, j]
        end
    end

    y .= A.y_full[A.fltni]
    return y
end

# Implement required LinearAlgebra interface
Base.size(A::StiffnessOperator) = (length(A.fltni), length(A.fltni))
Base.size(A::StiffnessOperator, d::Int) = size(A)[d]
Base.eltype(A::StiffnessOperator{T}) where {T} = T

"""
    stiffness_assembly(weights, H, Ht, dof_id) -> SparseMatrixCSC

Assemble global sparse stiffness matrix from metric weights.

Materializes K_el from metric weights (one-time cost), then assembles.
Used ONLY for building AMG preconditioner.
"""
function stiffness_assembly(
    weights::MetricWeightsAntiplane{T},
    H::Matrix{T},
    Ht::Matrix{T},
    dof_id::Array{Int,3}
) where {T}
    # Materialize K_el from metric weights (one-time cost for AMG setup)
    K_el = materialize_K_el_antiplane(weights, H, Ht)
    return stiffness_assembly(K_el, dof_id)
end

"""
    stiffness_assembly(K_el, dof_id) -> SparseMatrixCSC

Assemble global sparse stiffness matrix from elemental matrices.

Used ONLY for building AMG preconditioner (one-time cost).
"""
function stiffness_assembly(K_el::Array{T,3}, dof_id::Array{Int,3}) where {T}
    # Preallocate triplet vectors
    I = Vector{Int}(undef, length(K_el))
    J = Vector{Int}(undef, length(K_el))
    V = Vector{T}(undef, length(K_el))

    Nel = size(dof_id, 3)  # Number of elements
    ndof = maximum(dof_id)

    ct = 1
    for el in 1:Nel
        v = view(dof_id, :, :, el)
        for j in 1:length(v), i in 1:length(v)
            I[ct] = v[i]
            J[ct] = v[j]
            V[ct] = K_el[i, j, el]
            ct += 1
        end
    end

    # Construct sparse matrix (sum duplicates with '+' combiner)
    return sparse(I, J, V, ndof, ndof, +)
end


# ============================================================================
# Plane-strain matrix-free operator
# ============================================================================

"""
    apply_stiffness_plane_strain!(y, x, weights, H, Ht, dof_id, n_elements, ndof)

Apply plane-strain stiffness operator y = K*x using element-by-element tensor products.

The global vectors x, y have length 2*ndof in component-major order:
[u_x(1..ndof), u_y(1..ndof)].
"""
function apply_stiffness_plane_strain!(
    y::Vector{T},
    x::Vector{T},
    weights::MetricWeightsPlaneStrain{T},
    H::Matrix{T},
    Ht::Matrix{T},
    dof_id::Array{Int,3},
    n_elements::Int,
    ndof::Int
) where {T<:AbstractFloat}
    fill!(y, zero(T))

    nnodes = weights.nnodes
    ux_el = zeros(T, nnodes, nnodes)
    uy_el = zeros(T, nnodes, nnodes)
    fx_el = zeros(T, nnodes, nnodes)
    fy_el = zeros(T, nnodes, nnodes)
    tmp1  = zeros(T, nnodes, nnodes)
    tmp2  = zeros(T, nnodes, nnodes)

    for el in 1:n_elements
        local_idx = dof_id[:, :, el]

        # Gather: separate x and y components
        @inbounds for j in 1:nnodes, i in 1:nnodes
            idx = local_idx[i, j]
            ux_el[i, j] = x[idx]
            uy_el[i, j] = x[ndof + idx]
        end

        # Tensor-product element matvec
        apply_element_stiffness_plane_strain!(fx_el, fy_el, ux_el, uy_el,
                                              view(weights.g, :, :, :, el),
                                              H, Ht, tmp1, tmp2)

        # Scatter-add: both components
        @inbounds for j in 1:nnodes, i in 1:nnodes
            idx = local_idx[i, j]
            y[idx]        += fx_el[i, j]
            y[ndof + idx] += fy_el[i, j]
        end
    end

    return y
end

"""
    StiffnessOperatorPlaneStrain{T} <: AbstractMatrix{T}

Matrix-free stiffness operator for plane-strain formulation.

Operates on 2*ndof vectors in component-major order: [u_x(1..ndof), u_y(1..ndof)].
Uses tensor-product matvec with compact MetricWeightsPlaneStrain arrays.
"""
struct StiffnessOperatorPlaneStrain{T<:AbstractFloat} <: AbstractMatrix{T}
    weights::MetricWeightsPlaneStrain{T}
    H::Matrix{T}
    Ht::Matrix{T}
    dof_id::Array{Int,3}      # [nnodes, nnodes, n_elements] (spatial DOFs)
    n_elements::Int
    ndof::Int                  # spatial DOFs
    fltni::Vector{Int}         # free DOF indices in 2*ndof space
    x_full::Vector{T}          # workspace [2*ndof]
    y_full::Vector{T}          # workspace [2*ndof]
    # Element-local workspaces
    ux_el::Matrix{T}           # [nnodes, nnodes]
    uy_el::Matrix{T}           # [nnodes, nnodes]
    fx_el::Matrix{T}           # [nnodes, nnodes]
    fy_el::Matrix{T}           # [nnodes, nnodes]
    tmp1::Matrix{T}            # [nnodes, nnodes]
    tmp2::Matrix{T}            # [nnodes, nnodes]
end

function StiffnessOperatorPlaneStrain(
    weights::MetricWeightsPlaneStrain{T},
    H::Matrix{T},
    Ht::Matrix{T},
    dof_id::Array{Int,3},
    n_elements::Int,
    ndof::Int,
    fltni::Vector{Int}
) where {T<:AbstractFloat}
    x_full = zeros(T, 2 * ndof)
    y_full = zeros(T, 2 * ndof)
    nnodes = weights.nnodes
    ux_el = zeros(T, nnodes, nnodes)
    uy_el = zeros(T, nnodes, nnodes)
    fx_el = zeros(T, nnodes, nnodes)
    fy_el = zeros(T, nnodes, nnodes)
    tmp1  = zeros(T, nnodes, nnodes)
    tmp2  = zeros(T, nnodes, nnodes)
    return StiffnessOperatorPlaneStrain(weights, H, Ht, dof_id, n_elements, ndof, fltni,
                                        x_full, y_full,
                                        ux_el, uy_el, fx_el, fy_el, tmp1, tmp2)
end

function LinearAlgebra.mul!(y::Vector, A::StiffnessOperatorPlaneStrain{T}, x::Vector) where {T}
    ndof   = A.ndof
    nnodes = A.weights.nnodes

    fill!(A.x_full, zero(T))
    A.x_full[A.fltni] .= x

    fill!(A.y_full, zero(T))

    for el in 1:A.n_elements
        local_idx = A.dof_id[:, :, el]

        # Gather
        @inbounds for j in 1:nnodes, i in 1:nnodes
            idx = local_idx[i, j]
            A.ux_el[i, j] = A.x_full[idx]
            A.uy_el[i, j] = A.x_full[ndof + idx]
        end

        # Tensor-product matvec
        apply_element_stiffness_plane_strain!(A.fx_el, A.fy_el, A.ux_el, A.uy_el,
                                              view(A.weights.g, :, :, :, el),
                                              A.H, A.Ht, A.tmp1, A.tmp2)

        # Scatter-add
        @inbounds for j in 1:nnodes, i in 1:nnodes
            idx = local_idx[i, j]
            A.y_full[idx]        += A.fx_el[i, j]
            A.y_full[ndof + idx] += A.fy_el[i, j]
        end
    end

    y .= A.y_full[A.fltni]
    return y
end

Base.size(A::StiffnessOperatorPlaneStrain) = (length(A.fltni), length(A.fltni))
Base.size(A::StiffnessOperatorPlaneStrain, d::Int) = size(A)[d]
Base.eltype(A::StiffnessOperatorPlaneStrain{T}) where {T} = T

"""
    stiffness_assembly_plane_strain(weights, H, Ht, dof_id, ndof) -> SparseMatrixCSC

Assemble global sparse stiffness matrix for plane-strain (2*ndof × 2*ndof).

Materializes K_el from metric weights, then assembles sparse matrix.
Used ONLY for building the AMG preconditioner.
"""
function stiffness_assembly_plane_strain(
    weights::MetricWeightsPlaneStrain{T},
    H::Matrix{T},
    Ht::Matrix{T},
    dof_id::Array{Int,3},
    ndof::Int
) where {T}
    # Materialize K_el from metric weights (one-time cost for AMG setup)
    K_el = materialize_K_el_plane_strain(weights, H, Ht)
    return stiffness_assembly_plane_strain(K_el, dof_id, ndof)
end

"""
    stiffness_assembly_plane_strain(K_el, dof_id, ndof) -> SparseMatrixCSC

Assemble global sparse stiffness matrix for plane-strain (2*ndof × 2*ndof).

Used for building the AMG preconditioner.
"""
function stiffness_assembly_plane_strain(K_el::Array{T,3}, dof_id::Array{Int,3}, ndof::Int) where {T}
    Nel = size(dof_id, 3)
    nnodes_sq = size(dof_id, 1) * size(dof_id, 2)

    # Count total triplets: each element contributes (2*N)^2 entries
    total_entries = Nel * (2 * nnodes_sq)^2
    II = Vector{Int}(undef, total_entries)
    JJ = Vector{Int}(undef, total_entries)
    VV = Vector{T}(undef, total_entries)

    ct = 1
    for el in 1:Nel
        idx_spatial = dof_id[:, :, el][:]  # Spatial DOF indices, length N

        # Build global index mapping for this element: [idx_x; idx_y]
        global_idx = vcat(idx_spatial, ndof .+ idx_spatial)  # length 2N

        for j in 1:2*nnodes_sq, i in 1:2*nnodes_sq
            II[ct] = global_idx[i]
            JJ[ct] = global_idx[j]
            VV[ct] = K_el[i, j, el]
            ct += 1
        end
    end

    return sparse(II, JJ, VV, 2 * ndof, 2 * ndof, +)
end
