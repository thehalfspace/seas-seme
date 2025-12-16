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
    for i in 1:nnodes, j in 1:nnodes
        # Extract Jacobian components
        δx_δξ = jac_matrix[1, 1, i, j]  # ∂x/∂ξ
        δx_δη = jac_matrix[1, 2, i, j]  # ∂x/∂η
        δy_δξ = jac_matrix[2, 1, i, j]  # ∂y/∂ξ
        δy_δη = jac_matrix[2, 2, i, j]  # ∂y/∂η

        # Metric terms weighted by quadrature and material property
        # These represent the inverse metric tensor components
        denomξξ = δx_δξ * δx_δξ + δy_δξ * δy_δξ
        denomξη = δx_δξ * δx_δη + δy_δξ * δy_δη
        denomηξ = δx_δη * δx_δξ + δy_δη * δy_δξ
        denomηη = δx_δη * δx_δη + δy_δη * δy_δη

        # Check for division by zero and compute weights
        if denomξξ != 0
            Wξξ[i, j] = w2[i, j] * μ * jac_det[i, j] / denomξξ
        end
        if denomξη != 0
            Wξη[i, j] = w2[i, j] * μ * jac_det[i, j] / denomξη
        end
        if denomηξ != 0
            Wηξ[i, j] = w2[i, j] * μ * jac_det[i, j] / denomηξ
        end
        if denomηη != 0
            Wηη[i, j] = w2[i, j] * μ * jac_det[i, j] / denomηη
        end

        # Safety check for NaN/Inf
        if isnan(Wξξ[i, j]) || isinf(Wξξ[i, j])
            Wξξ[i, j] = 0.0
        end
        if isnan(Wξη[i, j]) || isinf(Wξη[i, j])
            Wξη[i, j] = 0.0
        end
        if isnan(Wηξ[i, j]) || isinf(Wηξ[i, j])
            Wηξ[i, j] = 0.0
        end
        if isnan(Wηη[i, j]) || isinf(Wηη[i, j])
            Wηη[i, j] = 0.0
        end
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

    return K_el
end
