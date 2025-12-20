"""
Stiffness matrix computation for spectral elements.

Implements elemental stiffness matrices for antiplane shear (SH waves).
"""

"""
    elemental_stiffness_matrix(μ::Real, jac_matrix::Array, jac_det::Matrix, basis) -> Matrix

Compute elemental stiffness matrix for antiplane shear in a spectral element.

# Arguments
- `μ::Real`: Shear modulus (Pa)
- `jac_matrix::Array`: Jacobian matrices [2, 2, nnodes, nnodes] for element
- `jac_det::Matrix`: Jacobian determinants [nnodes, nnodes] for element
- `basis`: Spectral element basis (LobattoLegendreBasis)

# Returns
- `Matrix{Float64}`: Elemental stiffness matrix [nnodes², nnodes²]

# Formulation
For antiplane shear (out-of-plane displacement w), the weak form involves:
```
K_e = ∫∫ (∂N/∂x ∂N/∂x + ∂N/∂y ∂N/∂y) μ dΩ
```

Transformed to reference coordinates (ξ, η):
```
K_e = ∑ w_p w_q [ (∂ξ/∂x)² + (∂ξ/∂y)² ] H_ip H_kp δ_jq δ_lq +
               [ (∂η/∂x)² + (∂η/∂y)² ] H_jq H_lq δ_ip δ_kp +
               cross terms
```

# Notes
- From Ampuero SEM Notes
- Includes NaN/Inf checks for numerical stability (division by zero)
- Matrix size is (nnodes² x nnodes²) to handle 2D element-wise indexing
"""
function elemental_stiffness_matrix(
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

    # Weight matrices with metric terms
    Wξξ = zeros(T, nnodes, nnodes)
    Wξη = zeros(T, nnodes, nnodes)
    Wηξ = zeros(T, nnodes, nnodes)
    Wηη = zeros(T, nnodes, nnodes)

    # Compute weight matrices incorporating Jacobian and shear modulus
    # Standard SEM form uses contravariant metric terms:
    #   G = J^{-1} J^{-T}
    # so weights are w_i w_j * μ * |J| * Gαβ
    for i in 1:nnodes, j in 1:nnodes
        # Jacobian entries: J = [x_ξ  x_η; y_ξ  y_η]
        x_ξ = jac_matrix[1, 1, i, j]
        x_η = jac_matrix[1, 2, i, j]
        y_ξ = jac_matrix[2, 1, i, j]
        y_η = jac_matrix[2, 2, i, j]

        detJ = jac_det[i, j]

        # Hard fail on invalid geometry; don't "fix" by clamping/zeroing
        if !isfinite(detJ) || abs(detJ) < 1e-18
            error("Bad Jacobian determinant at (i=$i, j=$j): detJ=$detJ")
        end
        if detJ < 0
            @warn "Negative Jacobian determinant (orientation flip?)" detJ i j
            # You can choose to use abs(detJ) instead if you want SPD no matter what:
            # detJ = abs(detJ)
        end

        inv_detJ = 1 / detJ

        # Inverse Jacobian entries: [ξ_x  ξ_y; η_x  η_y] = J^{-1}
        ξ_x =  y_η * inv_detJ
        ξ_y = -x_η * inv_detJ
        η_x = -y_ξ * inv_detJ
        η_y =  x_ξ * inv_detJ

        # Contravariant metric terms G = J^{-1} J^{-T}
        g11 = ξ_x*ξ_x + ξ_y*ξ_y
        g22 = η_x*η_x + η_y*η_y
        g12 = ξ_x*η_x + ξ_y*η_y

        wJμ = w2[i, j] * μ * detJ   # if detJ always > 0, this is fine
        # If you decide to force SPD even with negative detJ, use:
        # wJμ = w2[i, j] * μ * abs(detJ)

        Wξξ[i, j] = wJμ * g11
        Wηη[i, j] = wJμ * g22

        # Force symmetry: Wξη must equal Wηξ
        Wξη[i, j] = wJμ * g12
        Wηξ[i, j] = Wξη[i, j]
    end

    # Assemble stiffness matrix using tensor product structure
    Ke_temp = zeros(T, nnodes, nnodes, nnodes, nnodes)
    δ = Matrix{T}(I, nnodes, nnodes)  # Kronecker delta

    for i in 1:nnodes, j in 1:nnodes, k in 1:nnodes, l in 1:nnodes
        Ke_ξξ = zero(T)
        Ke_ξη = zero(T)
        Ke_ηξ = zero(T)
        Ke_ηη = zero(T)

        for p in 1:nnodes, q in 1:nnodes
            # Sum contributions from all quadrature points
            Ke_ξξ += Wξξ[p, q] * δ[j, q] * δ[l, q] * H[i, p] * H[k, p]
            Ke_ξη += Wξη[p, q] * δ[j, q] * δ[k, p] * H[i, p] * H[l, q]
            Ke_ηξ += Wηξ[p, q] * δ[i, p] * δ[l, q] * H[j, q] * H[k, p]
            Ke_ηη += Wηη[p, q] * δ[i, p] * δ[k, p] * H[j, q] * H[l, q]
        end

        # Sum all contributions
        Ke_temp[i, j, k, l] = Ke_ξξ + Ke_ηη + Ke_ξη + Ke_ηξ
    end

    # Reshape to 2D matrix (nnodes² x nnodes²)
    Ke = reshape(Ke_temp, nnodes * nnodes, nnodes * nnodes)

    return Ke
end

"""
    build_stiffness_matrices(mesh::UnstructuredSEMesh, material::MaterialProperties, basis)
        -> Array{Float64,3}

Build elemental stiffness matrices for all elements in the mesh.

# Arguments
- `mesh::UnstructuredSEMesh`: Complete spectral element mesh
- `material::MaterialProperties`: Material properties (shear modulus, etc.)
- `basis`: Spectral element basis

# Returns
- `K_el::Array{Float64,3}`: Elemental stiffness matrices [nnodes², nnodes², n_elements]

# Notes
- Stiffness matrices are NOT assembled globally (matrix-free approach)
- Each K_el[:,:,el] is a (nnodes² x nnodes²) matrix
- Used for both matrix-free operator and AMG preconditioner assembly

# Example
```julia
material = MaterialProperties(2670.0, 3464.0)
basis = LobattoLegendreBasis(4)
K_el = build_stiffness_matrices(mesh, material, basis)
```
"""
function build_stiffness_matrices(
    mesh::UnstructuredSEMesh{T},
    material::MaterialProperties{T},
    basis
) where {T<:AbstractFloat}
    nnodes = mesh.polynomial_degree + 1
    n_elements = mesh.n_elements

    # Allocate storage for elemental stiffness matrices
    K_el = zeros(T, nnodes * nnodes, nnodes * nnodes, n_elements)

    # Compute stiffness matrix for each element
    @info "Computing stiffness matrices for $(n_elements) elements..."
    for el in 1:n_elements
        K_el[:, :, el] = elemental_stiffness_matrix(
            material.μ,
            mesh.jac_matrix[:, :, :, :, el],
            mesh.jac_det[:, :, el],
            basis
        )
    end

    # Diagnostic: check for extreme values
    K_max = maximum(abs.(K_el))
    K_min = minimum(abs.(filter(!iszero, K_el)))
    @info "Stiffness matrix stats" K_max K_min material_μ=material.μ
    if K_max > 1e15 || K_min < 1e-10
        @warn "Extreme stiffness matrix values detected - likely bad Jacobians!"
        # Find which element has the problem
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
