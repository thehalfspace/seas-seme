"""
Element connectivity matrix construction for spectral element meshes.

This module builds the mapping from local (element-wise) node numbering
to global DOF numbering, ensuring C0 continuity across element boundaries
for Continuous Galerkin (CG) method.

Uses union-find (disjoint-set) data structure to correctly identify all
shared nodes across element interfaces, including corners and edges shared
by multiple elements.
"""

# ============================================================================
# Union-Find Data Structure
# ============================================================================

"""
    UnionFind

Disjoint-set data structure for tracking equivalence classes of local nodes.

Each local node is identified by a linear index computed from (element, i, j).
The union-find structure maintains parent pointers and compresses paths for
efficient queries.
"""
mutable struct UnionFind
    parent::Vector{Int}
    rank::Vector{Int}

    function UnionFind(n::Int)
        new(collect(1:n), zeros(Int, n))
    end
end

"""
    find!(uf::UnionFind, x::Int) -> Int

Find the root representative of element x with path compression.
"""
function find!(uf::UnionFind, x::Int)
    if uf.parent[x] != x
        uf.parent[x] = find!(uf, uf.parent[x])  # Path compression
    end
    return uf.parent[x]
end

"""
    union!(uf::UnionFind, x::Int, y::Int)

Union the sets containing x and y using union-by-rank.
"""
function union!(uf::UnionFind, x::Int, y::Int)
    root_x = find!(uf, x)
    root_y = find!(uf, y)

    if root_x == root_y
        return  # Already in same set
    end

    # Union by rank
    if uf.rank[root_x] < uf.rank[root_y]
        uf.parent[root_x] = root_y
    elseif uf.rank[root_x] > uf.rank[root_y]
        uf.parent[root_y] = root_x
    else
        uf.parent[root_y] = root_x
        uf.rank[root_x] += 1
    end
end

"""
    local_to_linear(el::Int, i::Int, j::Int, nnodes::Int, nel::Int) -> Int

Convert (element, i, j) local node coordinates to linear index.
"""
function local_to_linear(el::Int, i::Int, j::Int, nnodes::Int, nel::Int)
    return (el - 1) * nnodes * nnodes + (i - 1) * nnodes + j
end

"""
    linear_to_local(idx::Int, nnodes::Int, nel::Int) -> (Int, Int, Int)

Convert linear index back to (element, i, j) local coordinates.
"""
function linear_to_local(idx::Int, nnodes::Int, nel::Int)
    idx_0 = idx - 1
    el = div(idx_0, nnodes * nnodes) + 1
    remainder = idx_0 % (nnodes * nnodes)
    i = div(remainder, nnodes) + 1
    j = remainder % nnodes + 1
    return (el, i, j)
end

"""
    connectivity_matrix(mesh::UnstructuredMesh2D) -> Array{Int,3}

Build global DOF connectivity matrix for Continuous Galerkin SEM using union-find.

For each element `el`, `dof_id[:,:,el]` contains the global node IDs
corresponding to the local nodes in that element.

# Arguments
- `mesh::UnstructuredMesh2D`: Trixi.jl unstructured mesh

# Returns
- `Array{Int,3}`: Connectivity matrix `[nnodes, nnodes, n_elements]`

# Algorithm (Union-Find Approach)
1. Create union-find structure for all local nodes (nel × nnodes × nnodes)
2. For each interior interface in mesh.neighbour_information:
   - Extract face nodes from both elements (handling orientation)
   - Union corresponding node pairs across the interface
3. After all unions, enumerate unique components → assign global DOF IDs
4. Fill dof_id array with global IDs based on equivalence classes

# Notes
- Order-independent: processes all interfaces, not just reachable ones
- Correctly handles corners/edges shared by multiple elements
- Naturally parallelizable (union-find is thread-safe with atomic ops)
- Handles negative surface IDs (reversed orientation) via canonical_face_nodes
"""
function connectivity_matrix(mesh::UnstructuredMesh2D)::Array{Int,3}
    polydeg = mesh.polydeg
    nnodes = polydeg + 1
    nel = mesh.n_elements

    # Step 1: Create union-find structure for all local nodes
    total_local_nodes = nel * nnodes * nnodes
    uf = UnionFind(total_local_nodes)

    # Step 2: Process all interior interfaces and union shared nodes
    for i in 1:size(mesh.neighbour_information, 2)
        # Global edge orientation defined by corner nodes
        corner1_global = mesh.neighbour_information[1, i]  # Global corner ID (start)
        corner2_global = mesh.neighbour_information[2, i]  # Global corner ID (end)

        el1 = mesh.neighbour_information[3, i]
        el2 = mesh.neighbour_information[4, i]
        surf1_raw = mesh.neighbour_information[5, i]
        surf2_raw = mesh.neighbour_information[6, i]

        # Skip boundary faces (el2 == 0)
        if el2 == 0 || surf2_raw == 0
            continue
        end

        # Get face nodes from both elements in canonical order
        surf1_abs = abs(surf1_raw)
        surf2_abs = abs(surf2_raw)

        i1_idx, j1_idx = canonical_face_nodes(surf1_abs, nnodes)
        i2_idx, j2_idx = canonical_face_nodes(surf2_abs, nnodes)

        # Extract node indices for both faces in canonical order
        nodes1 = if i1_idx isa Int
            [(i1_idx, j) for j in j1_idx]
        else
            [(i, j1_idx) for i in i1_idx]
        end

        nodes2 = if i2_idx isa Int
            [(i2_idx, j) for j in j2_idx]
        else
            [(i, j2_idx) for i in i2_idx]
        end

        # Determine if canonical ordering matches global edge direction
        # by checking which corners are at the endpoints
        #
        # The mesh.element_node_ids gives us the corner nodes for each element
        el1_corners = mesh.element_node_ids[:, el1]  # 4 corners, counterclockwise
        el2_corners = mesh.element_node_ids[:, el2]

        # Map local face to its corner endpoints in canonical traversal order
        # Face 1 (bottom, η=-1): i varies 1→nnodes, corners are el_corners[1]→el_corners[2]
        # Face 2 (right, ξ=+1):  j varies 1→nnodes, corners are el_corners[2]→el_corners[3]
        # Face 3 (top, η=+1):    i varies 1→nnodes, corners are el_corners[4]→el_corners[3]
        # Face 4 (left, ξ=-1):   j varies 1→nnodes, corners are el_corners[1]→el_corners[4]
        local_corners1 = if surf1_abs == 1
            (el1_corners[1], el1_corners[2])
        elseif surf1_abs == 2
            (el1_corners[2], el1_corners[3])
        elseif surf1_abs == 3
            (el1_corners[4], el1_corners[3])
        else  # surf1_abs == 4
            (el1_corners[1], el1_corners[4])
        end

        local_corners2 = if surf2_abs == 1
            (el2_corners[1], el2_corners[2])
        elseif surf2_abs == 2
            (el2_corners[2], el2_corners[3])
        elseif surf2_abs == 3
            (el2_corners[4], el2_corners[3])
        else  # surf2_abs == 4
            (el2_corners[1], el2_corners[4])
        end

        # Determine if we need to reverse nodes to match global orientation
        # Global orientation is corner1_global → corner2_global
        #
        # Reverse if local traversal goes corner2→corner1 instead of corner1→corner2
        reverse1 = (local_corners1[1] == corner2_global && local_corners1[2] == corner1_global)
        reverse2 = (local_corners2[1] == corner2_global && local_corners2[2] == corner1_global)

        # Apply reversal if needed
        if reverse1
            nodes1 = reverse(nodes1)
        end
        if reverse2
            nodes2 = reverse(nodes2)
        end

        # Now both nodes1 and nodes2 are oriented consistently with the global edge
        # Union corresponding nodes
        for k in 1:nnodes
            i1, j1 = nodes1[k]
            i2, j2 = nodes2[k]

            linear1 = local_to_linear(el1, i1, j1, nnodes, nel)
            linear2 = local_to_linear(el2, i2, j2, nnodes, nel)

            union!(uf, linear1, linear2)
        end
    end

    # Step 3: Enumerate unique representatives and assign global DOF IDs
    representative_to_dof = Dict{Int, Int}()
    next_dof = 1

    dof_id = zeros(Int, nnodes, nnodes, nel)

    for linear_idx in 1:total_local_nodes
        root = find!(uf, linear_idx)

        if !haskey(representative_to_dof, root)
            representative_to_dof[root] = next_dof
            next_dof += 1
        end

        # Convert linear index back to (el, i, j) and assign DOF
        el, i, j = linear_to_local(linear_idx, nnodes, nel)
        dof_id[i, j, el] = representative_to_dof[root]
    end

    # Step 4: Validate CG connectivity
    total_dofs = maximum(dof_id)
    dg_dofs = nnodes^2 * nel  # Discontinuous Galerkin (no sharing)

    # For CG, expect 60-70% of DG DOFs for typical 2D quad meshes
    expected_cg_range = (0.6 * dg_dofs, 0.7 * dg_dofs)

    if total_dofs >= dg_dofs
        error("CG connectivity failed: total DOFs ($total_dofs) equals DG count ($dg_dofs). " *
              "Nodes are not being shared across element interfaces.")
    elseif total_dofs < expected_cg_range[1]
        @warn "DOF count suspiciously low" total_dofs expected_range=expected_cg_range
    elseif total_dofs > expected_cg_range[2]
        @warn "DOF count higher than expected for CG" total_dofs expected_range=expected_cg_range
    else
        @info "✓ CG connectivity successful: nodes shared across interfaces" total_dofs dg_dofs reduction=round((1 - total_dofs/dg_dofs)*100, digits=1)
    end

    print(dof_id)

    return dof_id
end

# Old sequential matching functions removed - replaced by union-find approach above
