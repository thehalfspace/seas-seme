"""
Boundary node extraction and impedance matrix computation.

This module identifies nodes on mesh boundaries (fault, creep, absorbing, free surface)
and computes their impedance matrix contributions for boundary conditions.
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
    identify_boundaries_geometric(mesh, boundary_name) -> boundary_el_id

Identify boundary elements by geometric location when boundary names are not preserved.

This function analyzes the corner coordinates of boundary elements to determine
which boundary they belong to based on spatial location.

# Arguments
- `mesh::UnstructuredMesh2D`: Trixi.jl mesh
- `boundary_name::Symbol`: Boundary to identify (:fault, :creep, :absorbing, :free_surface)

# Returns
- `boundary_el_id::Vector{CartesianIndex}`: Indices into mesh.neighbour_information

# Notes
Assumes standard fault geometry:
- Domain: x ∈ [0, Lx], y ∈ [0, Ly]
- Fault: right edge (x = Lx) from y_fault to Ly
- Creep: right edge (x = Lx) from 0 to y_fault
- Absorbing: left edge (x = 0) and bottom edge (y = 0)
- Free surface: top edge (y = Ly)
"""
function identify_boundaries_geometric(mesh::UnstructuredMesh2D, boundary_name::Symbol)
    # Get domain extents from corners
    x_coords = mesh.corners[1, :]
    y_coords = mesh.corners[2, :]

    x_min, x_max = extrema(x_coords)
    y_min, y_max = extrema(y_coords)

    # Tolerance for boundary identification (1% of domain size)
    tol_x = 0.01 * (x_max - x_min)
    tol_y = 0.01 * (y_max - y_min)

    # Estimate fault start (assume middle of right edge)
    y_fault = 0.5 * (y_min + y_max)

    boundary_el_id = CartesianIndex{2}[]

    # Loop through all boundary elements
    for i in 1:mesh.n_boundaries
        # Get the two corner nodes of this boundary edge
        node1_idx = mesh.neighbour_information[1, i]
        node2_idx = mesh.neighbour_information[2, i]

        x1, y1 = mesh.corners[1, node1_idx], mesh.corners[2, node1_idx]
        x2, y2 = mesh.corners[1, node2_idx], mesh.corners[2, node2_idx]

        # Midpoint of boundary edge
        x_mid = 0.5 * (x1 + x2)
        y_mid = 0.5 * (y1 + y2)

        # Determine which boundary this belongs to
        on_left = abs(x_mid - x_min) < tol_x
        on_right = abs(x_mid - x_max) < tol_x
        on_bottom = abs(y_mid - y_min) < tol_y
        on_top = abs(y_mid - y_max) < tol_y

        # Get element and surface IDs
        element_id = mesh.neighbour_information[4, i]
        surface_id = mesh.neighbour_information[5, i]

        # Classify boundary
        if boundary_name == :absorbing
            if on_left || on_bottom
                push!(boundary_el_id, CartesianIndex(surface_id, element_id))
            end

        elseif boundary_name == :free_surface
            if on_top
                push!(boundary_el_id, CartesianIndex(surface_id, element_id))
            end

        elseif boundary_name == :fault
            if on_right && y_mid >= y_fault - tol_y
                push!(boundary_el_id, CartesianIndex(surface_id, element_id))
            end

        elseif boundary_name == :creep
            if on_right && y_mid <= y_fault + tol_y
                push!(boundary_el_id, CartesianIndex(surface_id, element_id))
            end
        end
    end

    return boundary_el_id
end

"""
    get_boundary_nodes_structured(mesh, node_coords, jac_matrix, weights, impedance, dof_id, boundary_name)
        -> (node_ids, x_coords, y_coords, matrix)

Extract boundary nodes for structured meshes with axis-aligned boundaries.

# Arguments
- `mesh::UnstructuredMesh2D`: Trixi.jl mesh
- `node_coords::AbstractArray`: Node coordinates [2, nnodes, nnodes, nelements]
- `jac_matrix::AbstractArray`: Jacobian matrices [2, 2, nnodes, nnodes, nelements]
- `weights::AbstractVector`: Gauss-Lobatto quadrature weights
- `impedance::Real`: Impedance value (ρ*vs for absorbing, 1.0 for fault/creep)
- `dof_id::Array{Int,3}`: Connectivity matrix [nnodes, nnodes, nelements]
- `boundary_name::Symbol`: Boundary to extract (:fault, :creep, :absorbing, :free_surface)

# Returns
- `node_ids::Vector{Int}`: Unique global DOF IDs on boundary
- `x_coords::Vector`: x-coordinates of boundary nodes
- `y_coords::Vector`: y-coordinates of boundary nodes
- `matrix::Vector`: Impedance matrix contributions (summed for shared nodes)

# Notes
- For structured mesh, Right boundary is split: top half = fault, bottom half = creep
- Left and Bottom boundaries = absorbing
- Impedance matrix: weights[i] * jac1D[i] * impedance
- jac1D computed from Ampuero SEM notes (p. 23)
"""
function get_boundary_nodes_structured(
    mesh::UnstructuredMesh2D,
    node_coords::AbstractArray,
    jac_matrix::AbstractArray,
    weights::AbstractVector,
    impedance::Real,
    dof_id::Array{Int,3},
    boundary_name::Symbol
)
    polydeg = mesh.polydeg
    nnodes = polydeg + 1

    # Identify boundary elements based on name
    # For structured meshes: use symbolic boundary names
    # For unstructured meshes: use geometric location (boundary names may not be preserved)

    # Try to find boundaries by name first (works for structured meshes)
    if boundary_name == :absorbing
        # Try structured mesh names first
        boundary_el_id_left = findall(mesh.boundary_names .== :Left)
        boundary_el_id_bottom = findall(mesh.boundary_names .== :Bottom)
        boundary_el_id = unique(vcat(boundary_el_id_left, boundary_el_id_bottom))

        # If no named boundaries found, use geometric identification
        if isempty(boundary_el_id)
            boundary_el_id = identify_boundaries_geometric(mesh, :absorbing)
        end

    elseif boundary_name == :free_surface
        boundary_el_id = findall(mesh.boundary_names .== :Top)
        if isempty(boundary_el_id)
            boundary_el_id = identify_boundaries_geometric(mesh, :free_surface)
        end

    elseif boundary_name == :fault
        # For structured mesh: upper half of right boundary
        boundary_el_id_right = findall(mesh.boundary_names .== :Right)
        if !isempty(boundary_el_id_right)
            boundary_el_id = boundary_el_id_right[div(length(boundary_el_id_right), 2)+1:end]
        else
            boundary_el_id = identify_boundaries_geometric(mesh, :fault)
        end

    elseif boundary_name == :creep
        # For structured mesh: lower half of right boundary
        boundary_el_id_right = findall(mesh.boundary_names .== :Right)
        if !isempty(boundary_el_id_right)
            boundary_el_id = boundary_el_id_right[1:div(length(boundary_el_id_right), 2)]
        else
            boundary_el_id = identify_boundaries_geometric(mesh, :creep)
        end

    else
        error("Unknown boundary name: $boundary_name")
    end

    # Initialize storage
    boundary_node_id = Int[]
    boundary_x = Float64[]
    boundary_y = Float64[]
    boundary_mat = Float64[]
    boundary_mat_local = zeros(nnodes)

    # Loop over boundary elements
    for id in boundary_el_id
        surface = id[1]
        el = id[2]

        ncx = node_coords[1, :, :, el]
        ncy = node_coords[2, :, :, el]
        jac = jac_matrix[:, :, :, :, el]

        # Compute 1D Jacobian for boundary integration
        # From Ampuero SEM Notes page 23
        if surface == 1  # Bottom: sqrt((dx/dξ)² + (dy/dξ)²)
            jac1D = sqrt.(jac[1, 1, 1, :] .^ 2 .+ jac[2, 1, 1, :] .^ 2)

        elseif surface == 2  # Right: sqrt((dx/dη)² + (dy/dη)²)
            jac1D = sqrt.(jac[1, 2, :, end] .^ 2 .+ jac[2, 2, :, end] .^ 2)

        elseif surface == 3  # Top: sqrt((dx/dξ)² + (dy/dξ)²)
            jac1D = sqrt.(jac[1, 1, end, :] .^ 2 .+ jac[2, 1, end, :] .^ 2)

        elseif surface == 4  # Left: sqrt((dx/dη)² + (dy/dη)²)
            jac1D = sqrt.(jac[1, 2, :, 1] .^ 2 .+ jac[2, 2, :, 1] .^ 2)

        else
            error("Invalid surface ID: $surface")
        end

        # Impedance matrix contribution: weight * jac1D * impedance
        boundary_mat_local .= weights .* jac1D .* impedance

        # Extract boundary values
        append!(boundary_node_id, get_element_surface(dof_id[:, :, el], surface))
        append!(boundary_x, get_element_surface(ncx, surface))
        append!(boundary_y, get_element_surface(ncy, surface))
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
    get_element_surface(array::Array{Int,2}, surface::Int) -> Vector{Int}

Extract DOF IDs along a specified surface of an element.

# Arguments
- `array::Array{Int,2}`: Element DOF array [nnodes, nnodes]
- `surface::Int`: Surface ID (1=bottom, 2=right, 3=top, 4=left)

# Returns
- `Vector{Int}`: DOF IDs along surface

# Notes
- Surface numbering is counterclockwise from bottom
- Negative surface IDs indicate reversed orientation
"""
function get_element_surface(array::Array{Int,2}, surface::Int)::Vector{Int}
    if surface == 1  # Bottom
        return array[1, :]
    elseif surface == 2  # Right
        return array[:, end]
    elseif surface == 3  # Top
        return array[end, :]
    elseif surface == 4  # Left
        return array[:, 1]
    elseif surface == -1  # Bottom reversed
        return array[1, end:-1:1]
    elseif surface == -2  # Right reversed
        return array[end:-1:1, end]
    elseif surface == -3  # Top reversed
        return array[end, end:-1:1]
    elseif surface == -4  # Left reversed
        return array[end:-1:1, 1]
    else
        error("Invalid surface ID: $surface")
    end
end

"""
    get_element_surface(array::Array{T,2}, surface::Int) where T<:AbstractFloat -> Vector{T}

Extract coordinate values along a specified surface of an element.

# Arguments
- `array::Array{T,2}`: Element coordinate array [nnodes, nnodes]
- `surface::Int`: Surface ID (1=bottom, 2=right, 3=top, 4=left)

# Returns
- `Vector{T}`: Coordinate values along surface

# Notes
- Transposes array before extraction (coordinate convention)
- Surface numbering is counterclockwise from bottom
"""
function get_element_surface(array::Array{T,2}, surface::Int) where {T<:AbstractFloat}
    array_t = array'  # Transpose for coordinate convention

    if surface == 1  # Bottom
        return array_t[1, :]
    elseif surface == 2  # Right
        return array_t[:, end]
    elseif surface == 3  # Top
        return array_t[end, :]
    elseif surface == 4  # Left
        return array_t[:, 1]
    elseif surface == -1  # Bottom reversed
        return array_t[1, end:-1:1]
    elseif surface == -2  # Right reversed
        return array_t[end:-1:1, end]
    elseif surface == -3  # Top reversed
        return array_t[end, end:-1:1]
    elseif surface == -4  # Left reversed
        return array_t[end:-1:1, 1]
    else
        error("Invalid surface ID: $surface")
    end
end
