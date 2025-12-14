"""
Mass matrix computation for spectral elements.

Implements diagonal (lumped) mass matrix for explicit time integration.
"""

"""
    elemental_mass_matrix(ρ::Real, jac_det::Matrix, basis) -> Matrix

Compute elemental mass matrix for a single spectral element.

# Arguments
- `ρ::Real`: Material density (kg/m³)
- `jac_det::Matrix`: Jacobian determinants [nnodes x nnodes] for element
- `basis`: Spectral element basis (LobattoLegendreBasis from Trixi.jl)

# Returns
- `Matrix{Float64}`: Elemental mass matrix [nnodes x nnodes] (diagonal)

# Formula
```
M_e[i,j] = w[i] * w[j] * ρ * |J(ξᵢ, ηⱼ)|
```
where w are Gauss-Lobatto quadrature weights and J is the Jacobian.

# Notes
- Mass matrix is diagonal in spectral element method (mass lumping)
- Diagonal structure enables efficient explicit time stepping
- From Ampuero SEM Notes
"""
function elemental_mass_matrix(
    ρ::Real,
    jac_det::Matrix{T},
    basis
) where {T<:AbstractFloat}
    w = basis.weights
    Me = (w * w') .* ρ .* jac_det
    return Me
end

"""
    global_mass_matrix!(M_global, mesh, M_el, dof_id)

Assemble global mass matrix from elemental contributions.

# Arguments
- `M_global::Vector`: Global mass vector (modified in-place)
- `mesh::UnstructuredMesh2D`: Mesh
- `M_el::Array{T,3}`: Elemental mass matrices [nnodes, nnodes, n_elements]
- `dof_id::Array{Int,3}`: Connectivity matrix [nnodes, nnodes, n_elements]

# Returns
- `M_global::Vector`: Assembled global mass vector

# Notes
- Since mass matrix is diagonal, we store as a vector (not sparse matrix)
- Contributions from shared nodes (element boundaries) are summed
- Assembly loop is element-by-element
"""
function global_mass_matrix!(
    M_global::Vector{T},
    mesh::UnstructuredMesh2D,
    M_el::Array{T,3},
    dof_id::Array{Int,3}
) where {T<:AbstractFloat}
    for el in 1:mesh.n_elements
        local_index = dof_id[:, :, el]
        M_global[local_index] .+= M_el[:, :, el]
    end
    return M_global
end

"""
    build_mass_matrices(mesh::UnstructuredSEMesh, material::MaterialProperties, basis)
        -> (M_el, M_global)

Build both elemental and global mass matrices for the entire mesh.

# Arguments
- `mesh::UnstructuredSEMesh`: Complete spectral element mesh
- `material::MaterialProperties`: Material properties (density, etc.)
- `basis`: Spectral element basis

# Returns
- `M_el::Array{Float64,3}`: Elemental mass matrices [nnodes, nnodes, n_elements]
- `M_global::Vector{Float64}`: Global mass vector [ndof]

# Example
```julia
material = MaterialProperties(2670.0, 3464.0)  # ρ, vs
basis = LobattoLegendreBasis(4)
M_el, M_global = build_mass_matrices(mesh, material, basis)
```
"""
function build_mass_matrices(
    mesh::UnstructuredSEMesh{T},
    material::MaterialProperties{T},
    basis
) where {T<:AbstractFloat}
    nnodes = mesh.polynomial_degree + 1
    n_elements = mesh.n_elements

    # Allocate elemental mass matrices
    M_el = zeros(T, nnodes, nnodes, n_elements)

    # Compute elemental mass for each element
    for el in 1:n_elements
        M_el[:, :, el] = elemental_mass_matrix(
            material.ρ,
            mesh.jac_det[:, :, el],
            basis
        )
    end

    # Assemble global mass vector
    M_global = zeros(T, mesh.ndof)
    global_mass_matrix!(M_global, mesh.trixi_mesh, M_el, mesh.dof_id)

    return M_el, M_global
end
