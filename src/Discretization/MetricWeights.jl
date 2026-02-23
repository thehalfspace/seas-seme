"""
Metric weight arrays for tensor-product element stiffness computation.

Replaces stored K_el (full element stiffness matrices) with compact metric weight
arrays from which forces can be computed via tensor products. This is the standard
approach used by SPECFEM3D, NekRS, and sem2dpack.

Memory reduction:
- Old: K_el antiplane = [N², N², n_el] where N = nnodes² → O(N⁴) per element
- New: g antiplane = [nnodes, nnodes, 3, n_el] → O(N²) per element
- Reduction: ~40× for p=4 (nnodes=5), ~160× for plane-strain
"""

"""
    MetricWeightsAntiplane{T}

Compact metric weight arrays for antiplane shear (SH) elements.

Stores 3 scalar weight arrays per quadrature point:
  g[:,:,1,el] = Wξξ = w_i * w_j * μ * |J| * (ξ_x² + ξ_y²)
  g[:,:,2,el] = Wηη = w_i * w_j * μ * |J| * (η_x² + η_y²)
  g[:,:,3,el] = Wξη = w_i * w_j * μ * |J| * (ξ_x*η_x + ξ_y*η_y)  [= Wηξ by symmetry]

These match the 3 "elastic constants" arrays in the Fortran sem2dpack reference
(nelast=3 for SH general case: e1=Wξξ, e2=Wηη, e3=Wξη).

# Fields
- `g::Array{T,4}`: Metric weights [nnodes, nnodes, 3, n_elements]
- `n_elements::Int`: Number of elements
- `nnodes::Int`: Nodes per direction (polynomial_degree + 1)
"""
struct MetricWeightsAntiplane{T<:AbstractFloat, A<:AbstractArray{T,4}}
    g::A              # [nnodes, nnodes, 3, n_elements]
    n_elements::Int
    nnodes::Int
end

# Inner constructor for CPU path (default)
MetricWeightsAntiplane{T}(g::Array{T,4}, n_elements::Int, nnodes::Int) where {T<:AbstractFloat} =
    MetricWeightsAntiplane{T, Array{T,4}}(g, n_elements, nnodes)

"""
    MetricWeightsPlaneStrain{T}

Compact metric weight arrays for plane-strain P-SV elements.

Stores 10 scalar weight arrays per quadrature point, matching the Fortran sem2dpack
convention (nelast=10 for P-SV general case). The 10 components encode the
weighted contravariant metric terms multiplied by the constitutive coefficients
needed to compute K_xx, K_yy, and K_xy block contributions via tensor products.

Component layout (following sem2dpack a(:,:,1..10)):
  g[:,:,1,el]  = (λ+2μ) * wJ * ξ_x² + μ * wJ * ξ_y²    (K_xx, ξξ term)
  g[:,:,2,el]  = (λ+2μ) * wJ * η_x² + μ * wJ * η_y²    (K_xx, ηη term)
  g[:,:,3,el]  = (λ+2μ) * wJ * ξ_x*η_x + μ * wJ * ξ_y*η_y  (K_xx, ξη=ηξ term)
  g[:,:,4,el]  = μ * wJ * ξ_x² + (λ+2μ) * wJ * ξ_y²    (K_yy, ξξ term)
  g[:,:,5,el]  = μ * wJ * η_x² + (λ+2μ) * wJ * η_y²    (K_yy, ηη term)
  g[:,:,6,el]  = μ * wJ * ξ_x*η_x + (λ+2μ) * wJ * ξ_y*η_y  (K_yy, ξη=ηξ term)
  g[:,:,7,el]  = λ * wJ * ξ_x*ξ_y + μ * wJ * ξ_y*ξ_x   (K_xy, ξξ term)
  g[:,:,8,el]  = λ * wJ * η_x*η_y + μ * wJ * η_y*η_x   (K_xy, ηη term)
  g[:,:,9,el]  = λ * wJ * ξ_x*η_y + μ * wJ * ξ_y*η_x   (K_xy, ξη term)
  g[:,:,10,el] = λ * wJ * η_x*ξ_y + μ * wJ * η_y*ξ_x   (K_xy, ηξ term)

# Fields
- `g::Array{T,4}`: Metric weights [nnodes, nnodes, 10, n_elements]
- `n_elements::Int`: Number of elements
- `nnodes::Int`: Nodes per direction (polynomial_degree + 1)
"""
struct MetricWeightsPlaneStrain{T<:AbstractFloat, A<:AbstractArray{T,4}}
    g::A              # [nnodes, nnodes, 10, n_elements]
    n_elements::Int
    nnodes::Int
end

# Inner constructor for CPU path (default)
MetricWeightsPlaneStrain{T}(g::Array{T,4}, n_elements::Int, nnodes::Int) where {T<:AbstractFloat} =
    MetricWeightsPlaneStrain{T, Array{T,4}}(g, n_elements, nnodes)


"""
    compute_metric_weights_antiplane(μ, jac_matrix, jac_det, basis) -> (g, nnodes)

Compute metric weight arrays for one antiplane element.

Returns the 3 metric weight arrays Wξξ, Wηη, Wξη as a [nnodes, nnodes, 3] array.

# Arguments
- `μ::Real`: Shear modulus
- `jac_matrix::Array{T,4}`: Jacobian [2, 2, nnodes, nnodes]
- `jac_det::Matrix{T}`: Jacobian determinant [nnodes, nnodes]
- `basis`: LobattoLegendreBasis
"""
function compute_metric_weights_antiplane(
    μ::Real,
    jac_matrix::Array{T,4},
    jac_det::Matrix{T},
    basis
) where {T<:AbstractFloat}
    w = basis.weights
    w2 = w * w'  # Outer product of quadrature weights
    nnodes = length(w)

    g = zeros(T, nnodes, nnodes, 3)

    for i in 1:nnodes, j in 1:nnodes
        x_ξ = jac_matrix[1, 1, i, j]
        x_η = jac_matrix[1, 2, i, j]
        y_ξ = jac_matrix[2, 1, i, j]
        y_η = jac_matrix[2, 2, i, j]

        detJ = jac_det[i, j]

        if !isfinite(detJ) || abs(detJ) < 1e-18
            error("Bad Jacobian determinant at (i=$i, j=$j): detJ=$detJ")
        end
        if detJ < 0
            @warn "Negative Jacobian determinant (orientation flip?)" detJ i j
        end

        inv_detJ = 1 / detJ

        # Inverse Jacobian: [ξ_x ξ_y; η_x η_y] = J^{-1}
        ξ_x =  y_η * inv_detJ
        ξ_y = -x_η * inv_detJ
        η_x = -y_ξ * inv_detJ
        η_y =  x_ξ * inv_detJ

        # Contravariant metric G = J^{-1} J^{-T}
        g11 = ξ_x*ξ_x + ξ_y*ξ_y
        g22 = η_x*η_x + η_y*η_y
        g12 = ξ_x*η_x + ξ_y*η_y  # = g21 by symmetry

        wJμ = w2[i, j] * μ * detJ

        g[i, j, 1] = wJμ * g11  # Wξξ
        g[i, j, 2] = wJμ * g22  # Wηη
        g[i, j, 3] = wJμ * g12  # Wξη = Wηξ
    end

    return g
end


"""
    compute_metric_weights_plane_strain(λ, μ, jac_matrix, jac_det, basis) -> g

Compute metric weight arrays for one plane-strain element.

Returns the 10 metric weight arrays as a [nnodes, nnodes, 10] array.

# Arguments
- `λ::Real`: First Lamé parameter
- `μ::Real`: Shear modulus
- `jac_matrix::Array{T,4}`: Jacobian [2, 2, nnodes, nnodes]
- `jac_det::Matrix{T}`: Jacobian determinant [nnodes, nnodes]
- `basis`: LobattoLegendreBasis
"""
function compute_metric_weights_plane_strain(
    λ::Real,
    μ::Real,
    jac_matrix::Array{T,4},
    jac_det::Matrix{T},
    basis
) where {T<:AbstractFloat}
    w = basis.weights
    w2 = w * w'
    nnodes = length(w)

    g = zeros(T, nnodes, nnodes, 10)

    for i in 1:nnodes, j in 1:nnodes
        x_ξ = jac_matrix[1, 1, i, j]
        x_η = jac_matrix[1, 2, i, j]
        y_ξ = jac_matrix[2, 1, i, j]
        y_η = jac_matrix[2, 2, i, j]

        detJ = jac_det[i, j]

        if !isfinite(detJ) || abs(detJ) < 1e-18
            error("Bad Jacobian determinant at (i=$i, j=$j): detJ=$detJ")
        end
        if detJ < 0
            @warn "Negative Jacobian determinant (orientation flip?)" detJ i j
        end

        inv_detJ = 1 / detJ

        ξ_x =  y_η * inv_detJ
        ξ_y = -x_η * inv_detJ
        η_x = -y_ξ * inv_detJ
        η_y =  x_ξ * inv_detJ

        wJ = w2[i, j] * detJ
        λ2μ = λ + 2μ

        # K_xx components: (λ+2μ) ∂/∂x ∂/∂x + μ ∂/∂y ∂/∂y
        g[i, j, 1] = wJ * (λ2μ * ξ_x*ξ_x + μ * ξ_y*ξ_y)     # ξξ
        g[i, j, 2] = wJ * (λ2μ * η_x*η_x + μ * η_y*η_y)     # ηη
        g[i, j, 3] = wJ * (λ2μ * ξ_x*η_x + μ * ξ_y*η_y)     # ξη = ηξ (symmetric)

        # K_yy components: μ ∂/∂x ∂/∂x + (λ+2μ) ∂/∂y ∂/∂y
        g[i, j, 4] = wJ * (μ * ξ_x*ξ_x + λ2μ * ξ_y*ξ_y)     # ξξ
        g[i, j, 5] = wJ * (μ * η_x*η_x + λ2μ * η_y*η_y)     # ηη
        g[i, j, 6] = wJ * (μ * ξ_x*η_x + λ2μ * ξ_y*η_y)     # ξη = ηξ (symmetric)

        # K_xy components: λ ∂/∂x ∂/∂y + μ ∂/∂y ∂/∂x
        # Note: K_xy is NOT symmetric in the ξη/ηξ decomposition
        g[i, j, 7]  = wJ * (λ * ξ_x*ξ_y + μ * ξ_y*ξ_x)      # ξξ: both mixed-direction
        g[i, j, 8]  = wJ * (λ * η_x*η_y + μ * η_y*η_x)      # ηη: both mixed-direction
        g[i, j, 9]  = wJ * (λ * ξ_x*η_y + μ * ξ_y*η_x)      # ξη term
        g[i, j, 10] = wJ * (λ * η_x*ξ_y + μ * η_y*ξ_x)      # ηξ term

    end

    return g
end


"""
    build_metric_weights_antiplane(mesh, material, basis) -> MetricWeightsAntiplane

Build metric weight arrays for all antiplane elements.
"""
function build_metric_weights_antiplane(
    mesh::UnstructuredSEMesh{T},
    material::MaterialProperties{T},
    basis
) where {T<:AbstractFloat}
    nnodes = mesh.polynomial_degree + 1
    n_elements = mesh.n_elements

    g_all = zeros(T, nnodes, nnodes, 3, n_elements)

    @info "Computing antiplane metric weights for $(n_elements) elements..."
    for el in 1:n_elements
        g_all[:, :, :, el] = compute_metric_weights_antiplane(
            material.μ,
            mesh.jac_matrix[:, :, :, :, el],
            mesh.jac_det[:, :, el],
            basis
        )
    end

    # Diagnostic
    g_max = maximum(abs.(g_all))
    g_min = minimum(abs.(filter(!iszero, g_all)))
    @info "Antiplane metric weight stats" g_max g_min material_μ=material.μ
    if g_max > 1e15 || g_min < 1e-15
        @warn "Extreme metric weight values detected - likely bad Jacobians!"
    end

    return MetricWeightsAntiplane{T}(g_all, n_elements, nnodes)
end


"""
    build_metric_weights_plane_strain(mesh, material, basis) -> MetricWeightsPlaneStrain

Build metric weight arrays for all plane-strain elements.
"""
function build_metric_weights_plane_strain(
    mesh::UnstructuredSEMesh{T},
    material::MaterialProperties{T},
    basis
) where {T<:AbstractFloat}
    nnodes = mesh.polynomial_degree + 1
    n_elements = mesh.n_elements

    g_all = zeros(T, nnodes, nnodes, 10, n_elements)

    @info "Computing plane-strain metric weights for $(n_elements) elements..."
    for el in 1:n_elements
        g_all[:, :, :, el] = compute_metric_weights_plane_strain(
            material.λ,
            material.μ,
            mesh.jac_matrix[:, :, :, :, el],
            mesh.jac_det[:, :, el],
            basis
        )
    end

    # Diagnostic
    g_max = maximum(abs.(g_all))
    g_min = minimum(abs.(filter(!iszero, g_all)))
    @info "Plane-strain metric weight stats" g_max g_min material_λ=material.λ material_μ=material.μ
    if g_max > 1e15 || g_min < 1e-15
        @warn "Extreme metric weight values detected - likely bad Jacobians!"
    end

    return MetricWeightsPlaneStrain{T}(g_all, n_elements, nnodes)
end


"""
    apply_element_stiffness_antiplane!(f_el, u_el, g_el, H, Ht, tmp1, tmp2)

Compute f_el = K_el * u_el for one antiplane element using tensor products.

The tensor-product formulation (Ampuero SEM Notes, sem2dpack approach):
  f = H * (Wξξ .* (Hᵀ * u)) + (Wηη .* (u * H)) * Hᵀ
    + H * (Wξη .* (u * H)) + (Wξη .* (Hᵀ * u)) * Hᵀ

where H is the nnodes×nnodes derivative matrix (D_ij = dφ_j/dξ at node i),
u_el and f_el are nnodes×nnodes matrices (2D element DOF layout).

# Arguments
- `f_el::Matrix{T}`: Output force [nnodes, nnodes] (modified in-place)
- `u_el::Matrix{T}`: Input displacement [nnodes, nnodes]
- `g_el::Array{T,3}`: Metric weights [nnodes, nnodes, 3] for this element
- `H::Matrix{T}`: Derivative matrix [nnodes, nnodes]
- `Ht::Matrix{T}`: H' (pre-transposed)
- `tmp1::Matrix{T}`: Pre-allocated workspace [nnodes, nnodes]
- `tmp2::Matrix{T}`: Pre-allocated workspace [nnodes, nnodes]
"""
function apply_element_stiffness_antiplane!(
    f_el::Matrix{T},
    u_el::Matrix{T},
    g_el::AbstractArray{T,3},
    H::Matrix{T},
    Ht::Matrix{T},
    tmp1::Matrix{T},
    tmp2::Matrix{T}
) where {T<:AbstractFloat}
    # dU/dξ direction: Hᵀ * u_el  (nnodes×nnodes)
    mul!(tmp1, Ht, u_el)
    # dU/dη direction: u_el * H  (nnodes×nnodes)
    mul!(tmp2, u_el, H)

    # Wξξ contribution: H * (Wξξ .* dU/dξ)
    @inbounds for j in axes(tmp1, 2), i in axes(tmp1, 1)
        tmp1[i, j] *= g_el[i, j, 1]   # Wξξ * dU/dξ
    end

    # Wηη contribution: (Wηη .* dU/dη) * Hᵀ -> stored in f_el first
    @inbounds for j in axes(tmp2, 2), i in axes(tmp2, 1)
        tmp2[i, j] *= g_el[i, j, 2]   # Wηη * dU/dη
    end

    # f = H * tmp1 + tmp2 * Ht  (main diagonal-block contributions)
    mul!(f_el, H, tmp1)
    mul!(f_el, tmp2, Ht, one(T), one(T))  # f_el += tmp2 * Ht

    # Wξη cross-term 1: H * (Wξη .* dU/dη)
    # Need dU/dη again — recompute into tmp2 (was modified above)
    mul!(tmp2, u_el, H)   # restore dU/dη
    @inbounds for j in axes(tmp2, 2), i in axes(tmp2, 1)
        tmp2[i, j] *= g_el[i, j, 3]   # Wξη * dU/dη
    end
    mul!(f_el, H, tmp2, one(T), one(T))  # f_el += H * (Wξη * dU/dη)

    # Wξη cross-term 2: (Wξη .* dU/dξ) * Hᵀ
    # Need dU/dξ again — recompute into tmp1
    mul!(tmp1, Ht, u_el)  # restore dU/dξ
    @inbounds for j in axes(tmp1, 2), i in axes(tmp1, 1)
        tmp1[i, j] *= g_el[i, j, 3]   # Wξη * dU/dξ
    end
    mul!(f_el, tmp1, Ht, one(T), one(T))  # f_el += (Wξη * dU/dξ) * Hᵀ

    return f_el
end


"""
    apply_element_stiffness_plane_strain!(fx_el, fy_el, ux_el, uy_el, g_el, H, Ht, tmp1, tmp2)

Compute (fx_el, fy_el) = K_el * (ux_el, uy_el) for one plane-strain element.

Uses 10 metric weight components to evaluate the K_xx, K_yy, K_xy blocks via
tensor products. The result is:
  fx = K_xx * ux + K_xy * uy
  fy = K_xy' * ux + K_yy * uy

Each block is evaluated via the same tensor-product pattern as the antiplane case,
using the appropriate metric weight components.

# Arguments
- `fx_el, fy_el::Matrix{T}`: Output forces [nnodes, nnodes] (modified in-place)
- `ux_el, uy_el::Matrix{T}`: Input displacements [nnodes, nnodes]
- `g_el::Array{T,3}`: Metric weights [nnodes, nnodes, 10] for this element
- `H::Matrix{T}`: Derivative matrix [nnodes, nnodes]
- `Ht::Matrix{T}`: H' (pre-transposed)
- `tmp1, tmp2::Matrix{T}`: Pre-allocated workspaces [nnodes, nnodes]
"""
function apply_element_stiffness_plane_strain!(
    fx_el::Matrix{T},
    fy_el::Matrix{T},
    ux_el::Matrix{T},
    uy_el::Matrix{T},
    g_el::AbstractArray{T,3},
    H::Matrix{T},
    Ht::Matrix{T},
    tmp1::Matrix{T},
    tmp2::Matrix{T}
) where {T<:AbstractFloat}
    fill!(fx_el, zero(T))
    fill!(fy_el, zero(T))

    # ---- K_xx * ux -> contributes to fx ----
    # g components 1,2,3 (ξξ, ηη, ξη)
    _add_block_contribution!(fx_el, ux_el, g_el, 1, 2, 3, 3, H, Ht, tmp1, tmp2)

    # ---- K_yy * uy -> contributes to fy ----
    # g components 4,5,6 (ξξ, ηη, ξη)
    _add_block_contribution!(fy_el, uy_el, g_el, 4, 5, 6, 6, H, Ht, tmp1, tmp2)

    # ---- K_xy * uy -> contributes to fx ----
    # K_xy has asymmetric ξη/ηξ split → components 7,8,9,10
    _add_Kxy_contribution!(fx_el, uy_el, g_el, 7, 8, 9, 10, H, Ht, tmp1, tmp2)

    # ---- K_xy' * ux -> contributes to fy ----
    # K_yx = K_xy^T: swap roles of Hᵀ*u and u*H (transpose the operator)
    _add_Kxy_transpose_contribution!(fy_el, ux_el, g_el, 7, 8, 9, 10, H, Ht, tmp1, tmp2)

    return nothing
end


"""
    _add_block_contribution!(f, u, g_el, c_xi, c_eta, c_xieta, c_etaxi, H, Ht, tmp1, tmp2)

Add K_block * u contribution to f, where K_block has symmetric ξη/ηξ metric weights.

Used for K_xx (components 1,2,3,3) and K_yy (components 4,5,6,6).
"""
function _add_block_contribution!(
    f::Matrix{T}, u::Matrix{T}, g_el::AbstractArray{T,3},
    c_xi::Int, c_eta::Int, c_xieta::Int, c_etaxi::Int,
    H::Matrix{T}, Ht::Matrix{T}, tmp1::Matrix{T}, tmp2::Matrix{T}
) where T
    # dU/dξ = Hᵀ * u
    mul!(tmp1, Ht, u)
    # dU/dη = u * H
    mul!(tmp2, u, H)

    # Wξξ * dU/dξ -> H * result
    @inbounds for j in axes(tmp1, 2), i in axes(tmp1, 1)
        tmp1[i, j] *= g_el[i, j, c_xi]
    end
    mul!(f, H, tmp1, one(T), one(T))

    # Wηη * dU/dη -> result * Hᵀ
    @inbounds for j in axes(tmp2, 2), i in axes(tmp2, 1)
        tmp2[i, j] *= g_el[i, j, c_eta]
    end
    mul!(f, tmp2, Ht, one(T), one(T))

    # Wξη * dU/dη -> H * result  (cross-term 1)
    mul!(tmp2, u, H)   # restore dU/dη
    @inbounds for j in axes(tmp2, 2), i in axes(tmp2, 1)
        tmp2[i, j] *= g_el[i, j, c_xieta]
    end
    mul!(f, H, tmp2, one(T), one(T))

    # Wηξ * dU/dξ -> result * Hᵀ  (cross-term 2; for symmetric blocks c_etaxi == c_xieta)
    mul!(tmp1, Ht, u)  # restore dU/dξ
    @inbounds for j in axes(tmp1, 2), i in axes(tmp1, 1)
        tmp1[i, j] *= g_el[i, j, c_etaxi]
    end
    mul!(f, tmp1, Ht, one(T), one(T))

    return nothing
end


"""
    _add_Kxy_contribution!(f, u, g_el, c7, c8, c9, c10, H, Ht, tmp1, tmp2)

Add K_xy * u contribution to f.

K_xy = ∫ (λ ∂N/∂x ∂N/∂y + μ ∂N/∂y ∂N/∂x) dΩ
     = H * (g7 .* dU/dξ) + (g8 .* dU/dη) * Hᵀ
     + H * (g9 .* dU/dη) + (g10 .* dU/dξ) * Hᵀ

where g7..g10 are the 4 K_xy metric weight components.
"""
function _add_Kxy_contribution!(
    f::Matrix{T}, u::Matrix{T}, g_el::AbstractArray{T,3},
    c7::Int, c8::Int, c9::Int, c10::Int,
    H::Matrix{T}, Ht::Matrix{T}, tmp1::Matrix{T}, tmp2::Matrix{T}
) where T
    # dU/dξ = Hᵀ * u
    mul!(tmp1, Ht, u)
    # dU/dη = u * H
    mul!(tmp2, u, H)

    # g7 * dU/dξ -> H * result
    @inbounds for j in axes(tmp1, 2), i in axes(tmp1, 1)
        tmp1[i, j] *= g_el[i, j, c7]
    end
    mul!(f, H, tmp1, one(T), one(T))

    # g8 * dU/dη -> result * Hᵀ
    @inbounds for j in axes(tmp2, 2), i in axes(tmp2, 1)
        tmp2[i, j] *= g_el[i, j, c8]
    end
    mul!(f, tmp2, Ht, one(T), one(T))

    # g9 * dU/dη -> H * result
    mul!(tmp2, u, H)   # restore dU/dη
    @inbounds for j in axes(tmp2, 2), i in axes(tmp2, 1)
        tmp2[i, j] *= g_el[i, j, c9]
    end
    mul!(f, H, tmp2, one(T), one(T))

    # g10 * dU/dξ -> result * Hᵀ
    mul!(tmp1, Ht, u)  # restore dU/dξ
    @inbounds for j in axes(tmp1, 2), i in axes(tmp1, 1)
        tmp1[i, j] *= g_el[i, j, c10]
    end
    mul!(f, tmp1, Ht, one(T), one(T))

    return nothing
end


"""
    _add_Kxy_transpose_contribution!(f, u, g_el, c7, c8, c9, c10, H, Ht, tmp1, tmp2)

Add K_xy^T * u contribution to f.

K_xy^T = K_yx where K_yx[I,K] = K_xy[K,I].

Derivation shows that K_yx has the same tensor-product structure as K_xy in all
4 directional terms (ξξ, ηη, ξη, ηξ), with only the ξη and ηξ metric weight
components swapped: c9 ↔ c10.

This is because:
- ξξ term: g7 = wJ*(λ+μ)*ξ_x*ξ_y — same for K_xy and K_yx
- ηη term: g8 = wJ*(λ+μ)*η_x*η_y — same for K_xy and K_yx
- K_xy ξη uses g9 = wJ*(λ*ξ_x*η_y + μ*ξ_y*η_x), K_yx ξη uses g10
- K_xy ηξ uses g10 = wJ*(λ*η_x*ξ_y + μ*η_y*ξ_x), K_yx ηξ uses g9
"""
function _add_Kxy_transpose_contribution!(
    f::Matrix{T}, u::Matrix{T}, g_el::AbstractArray{T,3},
    c7::Int, c8::Int, c9::Int, c10::Int,
    H::Matrix{T}, Ht::Matrix{T}, tmp1::Matrix{T}, tmp2::Matrix{T}
) where T
    # K_yx uses the same pattern as K_xy but with c9 and c10 swapped
    _add_Kxy_contribution!(f, u, g_el, c7, c8, c10, c9, H, Ht, tmp1, tmp2)
    return nothing
end


"""
    materialize_K_el_antiplane(weights, H, Ht) -> Array{T,3}

Materialize full elemental stiffness matrices from metric weights.

Used ONLY for AMG preconditioner assembly (one-time cost at simulation start).
Returns K_el [nnodes², nnodes², n_elements] matching the old format.
"""
function materialize_K_el_antiplane(
    weights::MetricWeightsAntiplane{T},
    H::Matrix{T},
    Ht::Matrix{T}
) where T
    nnodes = weights.nnodes
    N = nnodes * nnodes
    n_elements = weights.n_elements

    K_el = zeros(T, N, N, n_elements)

    f_el  = zeros(T, nnodes, nnodes)
    tmp1  = zeros(T, nnodes, nnodes)
    tmp2  = zeros(T, nnodes, nnodes)
    u_el  = zeros(T, nnodes, nnodes)

    for el in 1:n_elements
        g_el = view(weights.g, :, :, :, el)
        for col in 1:N
            # Basis vector e_col (column of identity)
            fill!(u_el, zero(T))
            u_el[col] = one(T)

            fill!(f_el, zero(T))
            apply_element_stiffness_antiplane!(f_el, u_el, g_el, H, Ht, tmp1, tmp2)

            K_el[:, col, el] .= vec(f_el)
        end
    end

    return K_el
end


"""
    materialize_K_el_plane_strain(weights, H, Ht) -> Array{T,3}

Materialize full plane-strain elemental stiffness matrices from metric weights.

Returns K_el [2*nnodes², 2*nnodes², n_elements] matching the old format.
Used ONLY for AMG preconditioner assembly.
"""
function materialize_K_el_plane_strain(
    weights::MetricWeightsPlaneStrain{T},
    H::Matrix{T},
    Ht::Matrix{T}
) where T
    nnodes = weights.nnodes
    N = nnodes * nnodes
    n_elements = weights.n_elements

    K_el = zeros(T, 2N, 2N, n_elements)

    fx_el = zeros(T, nnodes, nnodes)
    fy_el = zeros(T, nnodes, nnodes)
    tmp1  = zeros(T, nnodes, nnodes)
    tmp2  = zeros(T, nnodes, nnodes)
    ux_el = zeros(T, nnodes, nnodes)
    uy_el = zeros(T, nnodes, nnodes)

    for el in 1:n_elements
        g_el = view(weights.g, :, :, :, el)

        # Columns 1..N correspond to ux DOFs
        for col in 1:N
            fill!(ux_el, zero(T)); ux_el[col] = one(T)
            fill!(uy_el, zero(T))
            apply_element_stiffness_plane_strain!(fx_el, fy_el, ux_el, uy_el, g_el, H, Ht, tmp1, tmp2)
            K_el[1:N,    col, el] .= vec(fx_el)
            K_el[N+1:2N, col, el] .= vec(fy_el)
        end

        # Columns N+1..2N correspond to uy DOFs
        for col in 1:N
            fill!(ux_el, zero(T))
            fill!(uy_el, zero(T)); uy_el[col] = one(T)
            apply_element_stiffness_plane_strain!(fx_el, fy_el, ux_el, uy_el, g_el, H, Ht, tmp1, tmp2)
            K_el[1:N,    N + col, el] .= vec(fx_el)
            K_el[N+1:2N, N + col, el] .= vec(fy_el)
        end
    end

    return K_el
end


# ============================================================================
# GPU helpers: upload/download metric weights
# ============================================================================

"""
    gpu(w::MetricWeightsAntiplane) -> MetricWeightsAntiplane with CuArray g

Upload metric weight array to GPU. Returns a new MetricWeightsAntiplane where
`g` is a CuArray. Requires CUDA.jl to be loaded.
"""
function gpu(w::MetricWeightsAntiplane{T}) where T
    g_cu = CuArray(w.g)
    return MetricWeightsAntiplane{T, typeof(g_cu)}(g_cu, w.n_elements, w.nnodes)
end

"""
    gpu(w::MetricWeightsPlaneStrain) -> MetricWeightsPlaneStrain with CuArray g
"""
function gpu(w::MetricWeightsPlaneStrain{T}) where T
    g_cu = CuArray(w.g)
    return MetricWeightsPlaneStrain{T, typeof(g_cu)}(g_cu, w.n_elements, w.nnodes)
end

"""
    cpu(w::MetricWeightsAntiplane) -> MetricWeightsAntiplane with CPU Array g

Download metric weight array from GPU to CPU.
"""
function cpu(w::MetricWeightsAntiplane{T}) where T
    g_cpu = Array(w.g)
    return MetricWeightsAntiplane{T, Array{T,4}}(g_cpu, w.n_elements, w.nnodes)
end

"""
    cpu(w::MetricWeightsPlaneStrain) -> MetricWeightsPlaneStrain with CPU Array g
"""
function cpu(w::MetricWeightsPlaneStrain{T}) where T
    g_cpu = Array(w.g)
    return MetricWeightsPlaneStrain{T, Array{T,4}}(g_cpu, w.n_elements, w.nnodes)
end
