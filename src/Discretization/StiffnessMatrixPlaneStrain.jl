"""
Plane-strain stiffness matrix computation for spectral elements.

Implements elemental stiffness matrices for in-plane deformation (P-SV waves)
with two displacement components (u_x, u_y).
"""

"""
    elemental_stiffness_matrix_plane_strain(λ, μ, jac_matrix, jac_det, basis) -> Matrix

Compute elemental stiffness matrix for plane-strain elasticity.

# Arguments
- `λ::Real`: First Lamé parameter (Pa)
- `μ::Real`: Shear modulus (Pa)
- `jac_matrix::Array`: Jacobian matrices [2, 2, nnodes, nnodes] for element
- `jac_det::Matrix`: Jacobian determinants [nnodes, nnodes] for element
- `basis`: Spectral element basis (LobattoLegendreBasis)

# Returns
- `Matrix{Float64}`: Elemental stiffness matrix [2*nnodes², 2*nnodes²]

# Formulation
For plane-strain with displacement (u_x, u_y), the constitutive matrix is:
```
C = [λ+2μ   λ     0]
    [λ     λ+2μ   0]
    [0      0      μ]
```

The element stiffness K_e = ∫∫ Bᵀ C B |J| dξ dη has block structure:
```
K_e = [K_xx  K_xy]
      [K_yx  K_yy]
```
where:
  K_xx[I,J] = ∫ ((λ+2μ) ∂Nᵢ/∂x ∂Nⱼ/∂x + μ ∂Nᵢ/∂y ∂Nⱼ/∂y) dΩ
  K_yy[I,J] = ∫ (μ ∂Nᵢ/∂x ∂Nⱼ/∂x + (λ+2μ) ∂Nᵢ/∂y ∂Nⱼ/∂y) dΩ
  K_xy[I,J] = ∫ (λ ∂Nᵢ/∂x ∂Nⱼ/∂y + μ ∂Nᵢ/∂y ∂Nⱼ/∂x) dΩ
  K_yx = K_xyᵀ

DOF ordering is component-major: [u_x(1..N), u_y(1..N)] where N = nnodes².

# Notes
- From Ampuero SEM Notes, extended to plane-strain
- Includes same Jacobian checks as antiplane stiffness
"""
function elemental_stiffness_matrix_plane_strain(
    λ::Real,
    μ::Real,
    jac_matrix::Array{T,4},
    jac_det::Matrix{T},
    basis
) where {T<:AbstractFloat}
    x = basis.nodes
    w = basis.weights
    H = basis.derivative_matrix'  # Transpose of derivative matrix

    w2 = w * w'  # Outer product of weights
    nnodes = length(x)
    N = nnodes * nnodes  # Spatial DOFs per element

    # Weight matrices for each block
    # We need separate metric-weighted matrices for each combination of
    # derivatives and constitutive coefficients.
    #
    # For K_xx: (λ+2μ) ∂/∂x ∂/∂x + μ ∂/∂y ∂/∂y
    # For K_yy: μ ∂/∂x ∂/∂x + (λ+2μ) ∂/∂y ∂/∂y
    # For K_xy: λ ∂/∂x ∂/∂y + μ ∂/∂y ∂/∂x
    #
    # In reference coordinates, ∂/∂x = ξ_x ∂/∂ξ + η_x ∂/∂η
    #                            ∂/∂y = ξ_y ∂/∂ξ + η_y ∂/∂η
    #
    # So each block has contributions from ξξ, ηη, ξη, and ηξ terms.

    # Pre-compute inverse Jacobian components and weighted metric terms
    # at each quadrature point
    ξ_x = zeros(T, nnodes, nnodes)
    ξ_y = zeros(T, nnodes, nnodes)
    η_x = zeros(T, nnodes, nnodes)
    η_y = zeros(T, nnodes, nnodes)
    wJ  = zeros(T, nnodes, nnodes)  # w_i * w_j * |J|

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

        ξ_x[i, j] =  y_η * inv_detJ
        ξ_y[i, j] = -x_η * inv_detJ
        η_x[i, j] = -y_ξ * inv_detJ
        η_y[i, j] =  x_ξ * inv_detJ
        wJ[i, j]  = w2[i, j] * detJ
    end

    # Allocate the four sub-blocks
    K_xx = zeros(T, N, N)
    K_yy = zeros(T, N, N)
    K_xy = zeros(T, N, N)

    # Build stiffness using tensor-product structure
    # Node (i,j) in element maps to local DOF index: I = (j-1)*nnodes + i
    # (column-major Julia ordering matching reshape convention)
    δ = Matrix{T}(I, nnodes, nnodes)

    # We assemble each block separately.
    # The key insight: shape function derivatives in reference space use
    # the same tensor-product structure as antiplane, but we have different
    # material coefficients for each combination.

    for i in 1:nnodes, j in 1:nnodes
        I = (j - 1) * nnodes + i
        for k in 1:nnodes, l in 1:nnodes
            K = (l - 1) * nnodes + k

            val_xx = zero(T)
            val_yy = zero(T)
            val_xy = zero(T)

            for p in 1:nnodes, q in 1:nnodes
                wJpq = wJ[p, q]

                # Physical derivatives of basis functions at (p,q):
                # ∂N_I/∂x = ξ_x * H[i,p]*δ[j,q] + η_x * δ[i,p]*H[j,q]
                # ∂N_I/∂y = ξ_y * H[i,p]*δ[j,q] + η_y * δ[i,p]*H[j,q]
                # Similarly for ∂N_K/∂x, ∂N_K/∂y

                dNI_dx = ξ_x[p, q] * H[i, p] * δ[j, q] + η_x[p, q] * δ[i, p] * H[j, q]
                dNI_dy = ξ_y[p, q] * H[i, p] * δ[j, q] + η_y[p, q] * δ[i, p] * H[j, q]
                dNK_dx = ξ_x[p, q] * H[k, p] * δ[l, q] + η_x[p, q] * δ[k, p] * H[l, q]
                dNK_dy = ξ_y[p, q] * H[k, p] * δ[l, q] + η_y[p, q] * δ[k, p] * H[l, q]

                # K_xx: (λ+2μ) ∂Ni/∂x ∂Nk/∂x + μ ∂Ni/∂y ∂Nk/∂y
                val_xx += wJpq * ((λ + 2μ) * dNI_dx * dNK_dx + μ * dNI_dy * dNK_dy)

                # K_yy: μ ∂Ni/∂x ∂Nk/∂x + (λ+2μ) ∂Ni/∂y ∂Nk/∂y
                val_yy += wJpq * (μ * dNI_dx * dNK_dx + (λ + 2μ) * dNI_dy * dNK_dy)

                # K_xy: λ ∂Ni/∂x ∂Nk/∂y + μ ∂Ni/∂y ∂Nk/∂x
                val_xy += wJpq * (λ * dNI_dx * dNK_dy + μ * dNI_dy * dNK_dx)
            end

            K_xx[I, K] = val_xx
            K_yy[I, K] = val_yy
            K_xy[I, K] = val_xy
        end
    end

    # Assemble full 2N×2N matrix in component-major order: [u_x; u_y]
    Ke = zeros(T, 2N, 2N)
    Ke[1:N, 1:N]         = K_xx
    Ke[1:N, N+1:2N]      = K_xy
    Ke[N+1:2N, 1:N]      = K_xy'   # K_yx = K_xy^T (symmetry)
    Ke[N+1:2N, N+1:2N]   = K_yy

    return Ke
end

"""
    build_stiffness_matrices_plane_strain(mesh, material, basis) -> Array{Float64,3}

Build elemental stiffness matrices for all elements (plane-strain formulation).

# Arguments
- `mesh::UnstructuredSEMesh`: Complete spectral element mesh
- `material::MaterialProperties`: Material properties (λ, μ)
- `basis`: Spectral element basis

# Returns
- `K_el::Array{Float64,3}`: Elemental stiffness matrices [2*nnodes², 2*nnodes², n_elements]
"""
function build_stiffness_matrices_plane_strain(
    mesh::UnstructuredSEMesh{T},
    material::MaterialProperties{T},
    basis
) where {T<:AbstractFloat}
    nnodes = mesh.polynomial_degree + 1
    n_elements = mesh.n_elements
    N = nnodes * nnodes  # Spatial DOFs per element

    # Allocate storage for elemental stiffness matrices
    K_el = zeros(T, 2N, 2N, n_elements)

    @info "Computing plane-strain stiffness matrices for $(n_elements) elements..."
    for el in 1:n_elements
        K_el[:, :, el] = elemental_stiffness_matrix_plane_strain(
            material.λ,
            material.μ,
            mesh.jac_matrix[:, :, :, :, el],
            mesh.jac_det[:, :, el],
            basis
        )
    end

    # Diagnostic: check for extreme values
    K_max = maximum(abs.(K_el))
    K_min = minimum(abs.(filter(!iszero, K_el)))
    @info "Plane-strain stiffness matrix stats" K_max K_min material_λ=material.λ material_μ=material.μ
    if K_max > 1e15 || K_min < 1e-10
        @warn "Extreme stiffness matrix values detected - likely bad Jacobians!"
        for el in 1:n_elements
            K_el_max = maximum(abs.(K_el[:, :, el]))
            if K_el_max > 1e15
                jac_det_el = mesh.jac_det[:, :, el]
                @warn "Element $el has extreme stiffness" K_el_max jac_det_min=minimum(jac_det_el) jac_det_max=maximum(jac_det_el)
                break
            end
        end
    end

    return K_el
end
