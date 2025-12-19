"""
Boundary node extraction and impedance matrix computation using Face abstraction.

This module provides robust boundary extraction that works correctly for
unstructured meshes with arbitrary element orientations.

# Key features
- Uses Face abstraction (element, local_face_id, flip) for orientation-agnostic geometry
- No manual Jacobian indexing or surface ID switching
- Geometric invariant checks available in debug mode
- Single code path for all boundary types
"""

"""
    BoundaryData{T}

Container for boundary node information.

# Fields
- `node_ids::Vector{Int}`: Global DOF indices of boundary nodes
- `coords::AbstractMatrix{T}`: Coordinates [2 x nnodes] (x and y)
- `matrix::Vector{T}`: Impedance matrix contributions for each node
"""
struct BoundaryData{T<:AbstractFloat}
    node_ids::Vector{Int}
    coords::AbstractMatrix{T}         # [2 x nnodes]
    matrix::Vector{T}
end

"""
    get_boundary_nodes(mesh, node_coords, jac_matrix, weights, impedance,
                      dof_id, boundary_name, nodes, face_map)
        -> (node_ids, x_coords, y_coords, matrix)

Extract boundary nodes using Face abstraction.

# Arguments
- `mesh::UnstructuredMesh2D`: Trixi.jl mesh
- `node_coords::AbstractArray`: Node coordinates [2, nnodes, nnodes, nelements]
- `jac_matrix::AbstractArray`: Jacobian matrices [2, 2, nnodes, nnodes, nelements]
- `weights::AbstractVector`: Gauss-Lobatto quadrature weights
- `impedance::Real`: Impedance value (ρ*vs for absorbing, 1.0 for fault/creep)
- `dof_id::Array{Int,3}`: Connectivity matrix [nnodes, nnodes, nelements]
- `boundary_name::Symbol`: Boundary to extract (:fault, :creep, :absorbing, :free_surface)
- `nodes::AbstractVector`: GLL node positions in [-1, 1]
- `face_map::FaceMap`: Prebuilt face map

# Returns
- `node_ids::Vector{Int}`: Unique global DOF IDs on boundary
- `x_coords::Vector`: x-coordinates of boundary nodes
- `y_coords::Vector`: y-coordinates of boundary nodes
- `matrix::Vector`: Impedance matrix contributions (summed for shared nodes)

# Algorithm
1. Get boundary faces from face_map
2. For each face, call face_geometry() to get J₁D
3. Compute impedance: weights[i] * J₁D[i] * impedance
4. Sum contributions for shared nodes

# Notes
- **No surface ID switching logic**
- **No manual Jacobian indexing**
- All geometry comes from face_geometry()
- Geometric checks enabled with SEAS_SEME_GEOMETRIC_CHECKS=1
"""
function get_boundary_nodes(
    mesh,
    node_coords::AbstractArray,
    jac_matrix::AbstractArray,
    weights::AbstractVector,
    impedance::Real,
    dof_id::Array{Int,3},
    boundary_name::Symbol,
    nodes::AbstractVector,
    face_map::FaceMap
)
    polydeg = mesh.polydeg
    nnodes = polydeg + 1

    # Get boundary face indices
    if !haskey(face_map.boundary_faces, boundary_name)
        error("Boundary '$boundary_name' not found in mesh. Available: $(keys(face_map.boundary_faces))")
    end

    boundary_face_indices = face_map.boundary_faces[boundary_name]

    # Initialize storage
    boundary_node_id = Int[]
    boundary_x = Float64[]
    boundary_y = Float64[]
    boundary_mat = Float64[]

    # Loop over boundary faces
    for idx in boundary_face_indices
        face_info = face_map.faces[idx]
        face = face_info.face
        el = face.element

        # Get element corners
        corners = mesh.corners[:, mesh.element_node_ids[:, el]]'

        # Get face geometry using the unbreakable abstraction
        geom = face_geometry(face, nodes, corners, jac_matrix[:, :, :, :, el])

        # Geometric check (if enabled)
        if get(ENV, "SEAS_SEME_GEOMETRIC_CHECKS", "0") == "1"
            try
                check_edge_length(geom, weights)
            catch e
                @warn "Edge length check failed for boundary $boundary_name, element $el, face $(face.local_face_id)" exception=e
            end
        end

        # Impedance matrix contribution: weight * J₁D * impedance
        boundary_mat_local = weights .* geom.J_1D .* impedance

        # Extract DOF IDs for this face
        face_node_ids = extract_face_values(dof_id[:, :, el], face)

        # Append to global lists
        append!(boundary_node_id, face_node_ids)
        append!(boundary_x, geom.x_phys)
        append!(boundary_y, geom.y_phys)
        append!(boundary_mat, boundary_mat_local)
    end

    # Sum contributions for shared nodes (edges/corners)
    boundary_node_id_unique = unique(boundary_node_id)
    boundary_mat_unique = Float64[]
    boundary_x_unique = Float64[]
    boundary_y_unique = Float64[]

    for n in boundary_node_id_unique
        id = findall(boundary_node_id .== n)
        push!(boundary_mat_unique, sum(boundary_mat[id]))
        push!(boundary_x_unique, boundary_x[id[1]])
        push!(boundary_y_unique, boundary_y[id[1]])
    end

    return boundary_node_id_unique, boundary_x_unique, boundary_y_unique, boundary_mat_unique
end
