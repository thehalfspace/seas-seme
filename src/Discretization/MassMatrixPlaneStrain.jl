"""
Plane-strain mass matrix computation for spectral elements.

For two-component displacement (u_x, u_y), the mass matrix is block-diagonal
with identical blocks for each component. The global mass vector has length
2*ndof in component-major order: [M_x(1..ndof), M_y(1..ndof)].
"""

"""
    build_mass_matrices_plane_strain(mesh, material, basis) -> (M_el, M_global)

Build elemental and global mass matrices for plane-strain formulation.

# Arguments
- `mesh::UnstructuredSEMesh`: Complete spectral element mesh
- `material::MaterialProperties`: Material properties (density)
- `basis`: Spectral element basis

# Returns
- `M_el::Array{Float64,3}`: Elemental mass matrices [nnodes, nnodes, n_elements]
  (same as antiplane — mass is per spatial node, identical for both components)
- `M_global::Vector{Float64}`: Global mass vector [2*ndof]
  Layout: [M_x(1..ndof), M_y(1..ndof)] where M_x == M_y

# Notes
The mass matrix for each component is the same scalar mass matrix as antiplane.
We store the duplicated form for simplicity so that vector operations on the
2*ndof state vectors work with elementwise division (a = f ./ M).
"""
function build_mass_matrices_plane_strain(
    mesh::UnstructuredSEMesh{T},
    material::MaterialProperties{T},
    basis
) where {T<:AbstractFloat}
    # Reuse the existing elemental mass matrix computation
    M_el, M_scalar = build_mass_matrices(mesh, material, basis)

    # Create 2*ndof global mass vector: [M_x; M_y]
    ndof = mesh.ndof
    M_global = zeros(T, 2 * ndof)
    M_global[1:ndof] .= M_scalar
    M_global[ndof+1:2*ndof] .= M_scalar

    return M_el, M_global
end
