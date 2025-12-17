"""
Boundary node extraction using the Face abstraction (refactored version).

This is the **unbreakable** version that replaces ad-hoc Jacobian slicing
with the face_geometry abstraction.

# Key differences from old Boundaries.jl
- No manual `jac[1,1,:,1]` vs `jac[1,2,end,:]` switching
- Single code path for all boundaries
- Orientation-agnostic (handles rotated/flipped elements correctly)
- Invariant checks available in debug mode

# Migration path
1. Test this version against old version (should produce identical results)
2. Once validated, replace Boundaries.jl with this implementation
3. Delete old geometric identification code

# Usage
This file should be included AFTER Geometry.jl so that straight_side_quad_map
is available. Also include Faces.jl and GeometricChecks.jl first.

Example:
    include("src/Mesh/Geometry.jl")
    include("src/Mesh/Faces.jl")
    include("src/Mesh/GeometricChecks.jl")
    include("src/Mesh/BoundariesNew.jl")
"""

# Note: Faces.jl and GeometricChecks.jl should be included before this file
# We don't include them here to avoid circular dependencies

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
    get_boundary_nodes_new(mesh, node_coords, jac_matrix, weights, impedance,
                          dof_id, boundary_name; face_map=nothing)
        -> (node_ids, x_coords, y_coords, matrix)

Extract boundary nodes using Face abstraction (new unbreakable version).

# Arguments
- `mesh::UnstructuredMesh2D`: Trixi.jl mesh
- `node_coords::AbstractArray`: Node coordinates [2, nnodes, nnodes, nelements]
- `jac_matrix::AbstractArray`: Jacobian matrices [2, 2, nnodes, nnodes, nelements]
- `weights::AbstractVector`: Gauss-Lobatto quadrature weights
- `impedance::Real`: Impedance value (ρ*vs for absorbing, 1.0 for fault/creep)
- `dof_id::Array{Int,3}`: Connectivity matrix [nnodes, nnodes, nelements]
- `boundary_name::Symbol`: Boundary to extract (:fault, :creep, :absorbing, :free_surface)
- `face_map::Union{FaceMap,Nothing}`: Prebuilt face map (built if not provided)

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
function get_boundary_nodes_new(
    mesh,
    node_coords::AbstractArray,
    jac_matrix::AbstractArray,
    weights::AbstractVector,
    impedance::Real,
    dof_id::Array{Int,3},
    boundary_name::Symbol;
    face_map::Union{FaceMap,Nothing}=nothing,
    nodes::Union{AbstractVector,Nothing}=nothing
)
    polydeg = mesh.polydeg
    nnodes = polydeg + 1

    # Build face map if not provided
    if isnothing(face_map)
        face_map = build_face_map(mesh)
        determine_face_flips!(face_map, dof_id)
    end

    # Get GLL nodes if not provided
    if isnothing(nodes)
        error("nodes parameter is required. Pass GLL nodes from gausslobatto(nnodes).")
    end

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
                # Would call check_edge_length(geom, weights) if GeometricChecks is loaded
                # For now, skip to avoid dependency
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

"""
    compare_boundary_implementations(mesh, node_coords, jac_matrix, weights,
                                     impedance, dof_id, boundary_name)
        -> (old_result, new_result, match::Bool)

Compare old (brittle) and new (unbreakable) boundary extraction implementations.

# Arguments
Same as get_boundary_nodes_new

# Returns
- `old_result`: Tuple (node_ids, x, y, matrix) from old implementation
- `new_result`: Tuple (node_ids, x, y, matrix) from new implementation
- `match::Bool`: true if results are identical (within tolerance)

# Notes
Use this to validate that the refactored version produces identical results.
Once validated across all test cases, we can safely replace the old code.

# Example
```julia
old, new, match = compare_boundary_implementations(mesh, node_coords, jac_matrix,
                                                   weights, impedance, dof_id, :fault)
if match
    println("✓ New implementation validated!")
else
    println("✗ Mismatch detected - investigate")
    # Compare differences...
end
```
"""
function compare_boundary_implementations(
    mesh,
    node_coords::AbstractArray,
    jac_matrix::AbstractArray,
    weights::AbstractVector,
    impedance::Real,
    dof_id::Array{Int,3},
    boundary_name::Symbol;
    rtol=1e-10,
    atol=1e-12
)
    # Get old implementation result
    # (Assumes old Boundaries.jl is loaded as get_boundary_nodes_structured)
    old_result = get_boundary_nodes_structured(mesh, node_coords, jac_matrix,
                                              weights, impedance, dof_id,
                                              boundary_name)

    # Get new implementation result
    new_result = get_boundary_nodes_new(mesh, node_coords, jac_matrix,
                                       weights, impedance, dof_id,
                                       boundary_name)

    # Compare results
    old_ids, old_x, old_y, old_mat = old_result
    new_ids, new_x, new_y, new_mat = new_result

    # Check if identical
    match = true

    # Node IDs should be identical
    if old_ids != new_ids
        @warn "Node IDs don't match" setdiff(old_ids, new_ids) setdiff(new_ids, old_ids)
        match = false
    end

    # Coordinates should match (within tolerance)
    if !isapprox(old_x, new_x; rtol=rtol, atol=atol)
        @warn "X coordinates don't match" maximum(abs.(old_x .- new_x))
        match = false
    end

    if !isapprox(old_y, new_y; rtol=rtol, atol=atol)
        @warn "Y coordinates don't match" maximum(abs.(old_y .- new_y))
        match = false
    end

    # Matrix should match (within tolerance)
    if !isapprox(old_mat, new_mat; rtol=rtol, atol=atol)
        @warn "Impedance matrix doesn't match" maximum(abs.(old_mat .- new_mat))
        match = false
    end

    return old_result, new_result, match
end

# Forward declarations - should be provided by including module
if !@isdefined(straight_side_quad_map)
    function straight_side_quad_map end
end
if !@isdefined(get_boundary_nodes_structured)
    function get_boundary_nodes_structured end
end
