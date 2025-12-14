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
- `coords::Matrix{T}`: Coordinates [2 x nnodes] (x and y)
- `matrix::Vector{T}`: Impedance matrix contributions for each node
"""
struct BoundaryData{T<:AbstractFloat}
    node_ids::Vector{Int}
    coords::Matrix{T}         # [2 x nnodes]
    matrix::Vector{T}
end

"""
    get_boundary_nodes_structured(mesh, node_coords, jac_matrix, weights, impedance, dof_id, boundary_name)
        -> (node_ids, x_coords, y_coords, matrix)

Extract boundary nodes for structured meshes with axis-aligned boundaries.

# Arguments
- `mesh::UnstructuredMesh2D`: Trixi.jl mesh
- `node_coords::Array{<:Any,4}`: Node coordinates [2, nnodes, nnodes, nelements]
- `jac_matrix::Array{<:Any,5}`: Jacobian matrices [2, 2, nnodes, nnodes, nelements]
- `weights::Vector`: Gauss-Lobatto quadrature weights
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
    weights::Vector,
    impedance::Real,
    dof_id::Array{Int,3},
    boundary_name::Symbol
)
    polydeg = mesh.polydeg
    nnodes = polydeg + 1

    # Identify boundary elements based on name
    if boundary_name == :absorbing
        # Left and Bottom boundaries are absorbing
        boundary_el_id = unique(
            vcat(
                findall(mesh.boundary_names .== :Left),
                findall(mesh.boundary_names .== :Bottom)
            )
        )

    elseif boundary_name == :free_surface
        boundary_el_id = findall(mesh.boundary_names .== :Top)

    elseif boundary_name == :fault
        # Right boundary, upper half (20-40 km depth)
        id_temp = findall(mesh.boundary_names .== :Right)
        boundary_el_id = id_temp[div(length(id_temp), 2)+1:end]

    elseif boundary_name == :creep
        # Right boundary, lower half (0-20 km depth)
        id_temp = findall(mesh.boundary_names .== :Right)
        boundary_el_id = id_temp[1:div(length(id_temp), 2)]

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
