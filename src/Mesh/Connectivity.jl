"""
Element connectivity matrix construction for spectral element meshes.

This module builds the mapping from local (element-wise) node numbering
to global DOF numbering, ensuring continuity across element boundaries.
"""

"""
    connectivity_matrix(mesh::UnstructuredMesh2D) -> Array{Int,3}

Build global DOF connectivity matrix from local element indices.

For each element `el`, `dof_id[:,:,el]` contains the global node IDs
corresponding to the local nodes in that element (organized as a 2D array
matching the spectral element grid).

# Arguments
- `mesh::UnstructuredMesh2D`: Trixi.jl unstructured mesh

# Returns
- `Array{Int,3}`: Connectivity matrix `[nnodes, nnodes, n_elements]`

# Algorithm
1. Initialize first element with sequential numbering 1:nnodes²
2. For each subsequent element:
   - Match boundary nodes with adjacent elements
   - Assign new global IDs to interior nodes
3. Ensure C0 continuity across interfaces

# Notes
- Surface numbering convention (counterclockwise from bottom):
  - Surface 1: Bottom edge (i=1, j=1:nnodes)
  - Surface 2: Right edge (i=1:nnodes, j=nnodes)
  - Surface 3: Top edge (i=nnodes, j=1:nnodes)
  - Surface 4: Left edge (i=1:nnodes, j=1)
- Negative surface IDs indicate reversed orientation
"""
function connectivity_matrix(mesh::UnstructuredMesh2D)::Array{Int,3}
    polydeg = mesh.polydeg
    nnodes = polydeg + 1
    nel = mesh.n_elements

    # Initialize DOF ID array
    dof_id = zeros(Int, nnodes, nnodes, nel)

    # First element: sequential numbering
    ig1 = reshape(collect(1:nnodes*nnodes), nnodes, nnodes)
    last_dof_id::Int = 0

    for el in 1:nel
        if el == 1
            dof_id[:, :, el] .= ig1
        end

        # Find surfaces (interfaces) for current element
        surf_id1 = findall(mesh.neighbour_information[3, :] .== el)
        surf_id2 = findall(mesh.neighbour_information[4, :] .== el)
        surf_id = sort(vcat(surf_id1, surf_id2))

        # Match boundaries with neighboring elements
        for surf in surf_id
            el1 = mesh.neighbour_information[3, surf]
            el2 = mesh.neighbour_information[4, surf]
            surf1 = mesh.neighbour_information[5, surf]
            surf2 = mesh.neighbour_information[6, surf]

            # Skip if this is a boundary surface (no neighbor)
            if el2 == 0 || surf2 == 0
                continue
            end

            # Ensure el1 is the element being processed (swap if needed)
            if any(dof_id[:, :, el1] .!= 0)
                el1, el2 = el2, el1
                surf1, surf2 = surf2, surf1
            end

            # Match DOF IDs based on surface orientations
            match_boundary_dofs!(dof_id, el1, el2, surf1, surf2, nnodes)
        end

        # Assign global IDs to remaining (interior) nodes
        indx = findall(dof_id[:, :, el] .== 0)
        for i in indx
            dof_id[i, el] = last_dof_id + 1
            last_dof_id += 1
        end

        last_dof_id = dof_id[nnodes, nnodes, el]
    end

    return dof_id
end

"""
    match_boundary_dofs!(dof_id, el1, el2, surf1, surf2, nnodes)

Match boundary DOFs between two adjacent elements based on surface orientations.

# Arguments
- `dof_id::Array{Int,3}`: Connectivity matrix (modified in-place)
- `el1::Int`: Element 1 ID (to be filled)
- `el2::Int`: Element 2 ID (already filled)
- `surf1::Int`: Surface ID on element 1
- `surf2::Int`: Surface ID on element 2
- `nnodes::Int`: Number of nodes per element edge

# Notes
Surface numbering:
- 1: Bottom (i=1)
- 2: Right (j=nnodes)
- 3: Top (i=nnodes)
- 4: Left (j=1)
Negative values indicate reversed orientation.
"""
function match_boundary_dofs!(dof_id::Array{Int,3}, el1::Int, el2::Int,
                               surf1::Int, surf2::Int, nnodes::Int)
    # Surface 1 (bottom) of el1
    if surf1 == 1
        if surf2 == 1
            dof_id[1, :, el1] = dof_id[1, :, el2]
        elseif surf2 == 2
            dof_id[1, :, el1] = dof_id[:, nnodes, el2]
        elseif surf2 == 3
            dof_id[1, :, el1] = dof_id[nnodes, :, el2]
        elseif surf2 == 4
            dof_id[1, :, el1] = dof_id[:, 1, el2]
        elseif surf2 == -1
            dof_id[1, :, el1] = dof_id[1, end:-1:1, el2]
        elseif surf2 == -2
            dof_id[1, :, el1] = dof_id[end:-1:1, nnodes, el2]
        elseif surf2 == -3
            dof_id[1, :, el1] = dof_id[nnodes, end:-1:1, el2]
        elseif surf2 == -4
            dof_id[1, :, el1] = dof_id[end:-1:1, 1, el2]
        end

    # Surface 2 (right) of el1
    elseif surf1 == 2
        if surf2 == 1
            dof_id[:, nnodes, el1] = dof_id[1, :, el2]
        elseif surf2 == 2
            dof_id[:, nnodes, el1] = dof_id[:, nnodes, el2]
        elseif surf2 == 3
            dof_id[:, nnodes, el1] = dof_id[nnodes, :, el2]
        elseif surf2 == 4
            dof_id[:, nnodes, el1] = dof_id[:, 1, el2]
        elseif surf2 == -1
            dof_id[:, nnodes, el1] = dof_id[1, end:-1:1, el2]
        elseif surf2 == -2
            dof_id[:, nnodes, el1] = dof_id[end:-1:1, nnodes, el2]
        elseif surf2 == -3
            dof_id[:, nnodes, el1] = dof_id[nnodes, end:-1:1, el2]
        elseif surf2 == -4
            dof_id[:, nnodes, el1] = dof_id[end:-1:1, 1, el2]
        end

    # Surface 3 (top) of el1
    elseif surf1 == 3
        if surf2 == 1
            dof_id[nnodes, :, el1] = dof_id[1, :, el2]
        elseif surf2 == 2
            dof_id[nnodes, :, el1] = dof_id[:, nnodes, el2]
        elseif surf2 == 3
            dof_id[nnodes, :, el1] = dof_id[nnodes, :, el2]
        elseif surf2 == 4
            dof_id[nnodes, :, el1] = dof_id[:, 1, el2]
        elseif surf2 == -1
            dof_id[nnodes, :, el1] = dof_id[1, end:-1:1, el2]
        elseif surf2 == -2
            dof_id[nnodes, :, el1] = dof_id[end:-1:1, nnodes, el2]
        elseif surf2 == -3
            dof_id[nnodes, :, el1] = dof_id[nnodes, end:-1:1, el2]
        elseif surf2 == -4
            dof_id[nnodes, :, el1] = dof_id[end:-1:1, 1, el2]
        end

    # Surface 4 (left) of el1
    elseif surf1 == 4
        if surf2 == 1
            dof_id[:, 1, el1] = dof_id[1, :, el2]
        elseif surf2 == 2
            dof_id[:, 1, el1] = dof_id[:, nnodes, el2]
        elseif surf2 == 3
            dof_id[:, 1, el1] = dof_id[nnodes, :, el2]
        elseif surf2 == 4
            dof_id[:, 1, el1] = dof_id[:, 1, el2]
        elseif surf2 == -1
            dof_id[:, 1, el1] = dof_id[1, end:-1:1, el2]
        elseif surf2 == -2
            dof_id[:, 1, el1] = dof_id[end:-1:1, nnodes, el2]
        elseif surf2 == -3
            dof_id[:, 1, el1] = dof_id[nnodes, end:-1:1, el2]
        elseif surf2 == -4
            dof_id[:, 1, el1] = dof_id[end:-1:1, 1, el2]
        end
    end

    return nothing
end
