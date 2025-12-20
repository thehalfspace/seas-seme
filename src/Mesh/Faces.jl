"""
Face abstraction for spectral element boundary operations.

This module provides a single source of truth for face geometry, eliminating
ad-hoc indexing into Jacobian arrays and orientation guesswork.

# Design Philosophy

Instead of "guessing which reference side you're on" (surface=1..4, slices of jac[...]),
we build a **single source of truth** for every face:

    Face = (element id, local face id, orientation)

Everything (nodes, tangents, normals, J₁D, lifting/penalty ops) is derived from that one record.

# References
- Battle-tested patterns from Nektar++, deal.II, MFEM, Trixi-style DG
- Ampuero SEM Notes (page 23) for boundary Jacobian formulas
"""

"""
    Face

Single source of truth for a mesh face.

# Fields
- `element::Int`: Element ID (1-based)
- `local_face_id::Int`: Canonical local face ID ∈ {1,2,3,4}
  - 1: η = -1 (bottom), nodes vary in ξ
  - 2: ξ = +1 (right),  nodes vary in η
  - 3: η = +1 (top),    nodes vary in ξ
  - 4: ξ = -1 (left),   nodes vary in η
- `flip::Bool`: Orientation flag (true = reversed node ordering)

# Notes
Unstructured quads are the same reference element rotated/flipped.
If we encode orientation once per face, we can reuse the same extraction logic safely.
"""
struct Face
    element::Int
    local_face_id::Int  # 1-4 canonical
    flip::Bool          # orientation relative to canonical
end

"""
    FaceType

Classification of face location in the mesh.

# Values
- `Interior`: Face connects two elements
- `Boundary`: Face on domain boundary
"""
@enum FaceType Interior Boundary

"""
    FaceInfo

Complete information about a face including its neighbors and physical properties.

# Fields
- `face::Face`: Primary face descriptor
- `type::FaceType`: Interior or Boundary
- `neighbor_element::Int`: Neighbor element ID (0 for boundary)
- `neighbor_local_face_id::Int`: Neighbor's local face ID (0 for boundary)
- `neighbor_flip::Bool`: Neighbor's orientation flag
- `boundary_name::Symbol`: Boundary name if type==Boundary (e.g., :fault, :absorbing)

# Notes
This is the complete record for a face. Most operations only need the `face` field.
"""
struct FaceInfo
    face::Face
    type::FaceType
    neighbor_element::Int
    neighbor_local_face_id::Int
    neighbor_flip::Bool
    boundary_name::Symbol

    # Constructor for boundary faces
    function FaceInfo(face::Face, boundary_name::Symbol)
        new(face, Boundary, 0, 0, false, boundary_name)
    end

    # Constructor for interior faces
    function FaceInfo(face::Face, neighbor_element::Int,
                     neighbor_local_face_id::Int, neighbor_flip::Bool)
        new(face, Interior, neighbor_element, neighbor_local_face_id,
            neighbor_flip, :none)
    end
end

"""
    FaceMap

Container for all face information in the mesh.

# Fields
- `faces::Vector{FaceInfo}`: All faces (interior + boundary)
- `boundary_faces::Dict{Symbol, Vector{Int}}`: Indices into `faces` by boundary name

# Notes
Build once at mesh setup, then query for all boundary operations.
"""
struct FaceMap
    faces::Vector{FaceInfo}
    boundary_faces::Dict{Symbol, Vector{Int}}
end

"""
    canonical_face_nodes(local_face_id::Int, nnodes::Int) -> (i_indices, j_indices)

Get reference element node indices for a canonical face.

# Arguments
- `local_face_id::Int`: Face ID ∈ {1,2,3,4}
- `nnodes::Int`: Number of nodes per element edge (polydeg + 1)

# Returns
- `(i_indices, j_indices)`: Tuple of index ranges for ξ and η directions

# Face conventions
- Face 1 (bottom): η = -1, j = 1,      i = 1:nnodes
- Face 2 (right):  ξ = +1, i = nnodes, j = 1:nnodes
- Face 3 (top):    η = +1, j = nnodes, i = 1:nnodes
- Face 4 (left):   ξ = -1, i = 1,      j = 1:nnodes

# Example
```julia
i_idx, j_idx = canonical_face_nodes(2, 5)  # Right face, polydeg=4
# Returns: (5, 1:5) -> nodes at ξ=+1, varying η
```
"""
function canonical_face_nodes(local_face_id::Int, nnodes::Int)
    if local_face_id == 1  # Bottom: η = -1
        return (1:nnodes, 1)
    elseif local_face_id == 2  # Right: ξ = +1
        return (nnodes, 1:nnodes)
    elseif local_face_id == 3  # Top: η = +1
        return (1:nnodes, nnodes)
    elseif local_face_id == 4  # Left: ξ = -1
        return (1, 1:nnodes)
    else
        error("Invalid local_face_id: $local_face_id (must be 1-4)")
    end
end

"""
    reference_coordinates(local_face_id::Int, nodes::AbstractVector) -> (ξ_vals, η_vals)

Get reference coordinates (ξ, η) for nodes along a face.

# Arguments
- `local_face_id::Int`: Face ID ∈ {1,2,3,4}
- `nodes::AbstractVector`: GLL node positions in [-1, 1]

# Returns
- `(ξ_vals, η_vals)`: Vectors of reference coordinates for each face node

# Notes
Returns coordinates in canonical order. Apply `reverse` if face.flip == true.

# Example
```julia
nodes = [-1.0, -0.5, 0.0, 0.5, 1.0]  # 5 GLL points
ξ, η = reference_coordinates(2, nodes)  # Right face
# ξ = [1.0, 1.0, 1.0, 1.0, 1.0]  (constant at ξ=+1)
# η = [-1.0, -0.5, 0.0, 0.5, 1.0]  (varying)
```
"""
function reference_coordinates(local_face_id::Int, nodes::AbstractVector)
    nnodes = length(nodes)

    if local_face_id == 1  # Bottom: η = -1
        return (nodes, fill(-1.0, nnodes))
    elseif local_face_id == 2  # Right: ξ = +1
        return (fill(1.0, nnodes), nodes)
    elseif local_face_id == 3  # Top: η = +1
        return (nodes, fill(1.0, nnodes))
    elseif local_face_id == 4  # Left: ξ = -1
        return (fill(-1.0, nnodes), nodes)
    else
        error("Invalid local_face_id: $local_face_id")
    end
end

"""
    extract_face_values(array::AbstractMatrix, face::Face) -> Vector

Extract values along a face from a 2D element array.

# Arguments
- `array::AbstractMatrix`: Element data [nnodes, nnodes]
- `face::Face`: Face descriptor

# Returns
- `Vector`: Values along the face (reversed if face.flip == true)

# Notes
This replaces all the ad-hoc `array[1,:]`, `array[:,end]` indexing.
Works for both DOF IDs and coordinate arrays.

# Example
```julia
dof_array = reshape(1:25, 5, 5)  # 5x5 element
face = Face(1, 2, false)  # Right face, no flip
ids = extract_face_values(dof_array, face)
# Returns: [5, 10, 15, 20, 25] (right column)
```
"""
function extract_face_values(array::AbstractMatrix, face::Face)
    i_idx, j_idx = canonical_face_nodes(face.local_face_id, size(array, 1))

    # Extract values
    if i_idx isa Int
        values = array[i_idx, j_idx]
    else
        values = array[i_idx, j_idx]
    end

    # Apply flip if needed
    return face.flip ? reverse(values) : values
end

# Special method for coordinate arrays (requires transpose)
function extract_face_values(array::AbstractMatrix{T}, face::Face) where {T<:AbstractFloat}
    array_t = array'  # Coordinate convention
    i_idx, j_idx = canonical_face_nodes(face.local_face_id, size(array_t, 1))

    # Extract values
    if i_idx isa Int
        values = array_t[i_idx, j_idx]
    else
        values = array_t[i_idx, j_idx]
    end

    # Apply flip if needed
    return face.flip ? reverse(values) : values
end

"""
    face_geometry(face::Face, nodes::AbstractVector, corners::AbstractMatrix,
                  jac_matrix::AbstractArray)
        -> (x_phys, y_phys, t_x, t_y, n_x, n_y, J_1D)

Compute complete geometric information for a face.

# Arguments
- `face::Face`: Face descriptor
- `nodes::AbstractVector`: GLL node positions in [-1, 1]
- `corners::AbstractMatrix`: Element corner points [4 x 2]
- `jac_matrix::AbstractArray`: Jacobian matrix [2, 2, nnodes, nnodes]

# Returns
Named tuple with:
- `x_phys::Vector`: Physical x-coordinates of face nodes
- `y_phys::Vector`: Physical y-coordinates of face nodes
- `t_x::Vector`: Tangent vector x-components
- `t_y::Vector`: Tangent vector y-components
- `n_x::Vector`: Outward normal x-components (NOT normalized)
- `n_y::Vector`: Outward normal y-components (NOT normalized)
- `J_1D::Vector`: Line Jacobian (for integration: ∫ f ds = ∑ w_i J_1D[i] f[i])

# Algorithm
For each face quadrature point:
1. Evaluate physical coordinates x_f(s) from element mapping
2. Evaluate tangent t = dx/ds from Jacobian matrix
3. Compute line Jacobian J₁D = |t|
4. Compute outward normal n = (t_y, -t_x) / |t| (scaled by |t| for unnormalized version)

# Notes
- **Orientation-agnostic**: Doesn't care if face is "ξ-constant" or "η-constant"
- Computes what the edge **is** in physical space
- Normals point outward from the element
- Flip is handled automatically via reference_coordinates
- This replaces all the `if surface == 1 ... elseif surface == 2 ...` logic

# Invariant check
The following should hold (within tolerance):
```julia
geom = face_geometry(...)
edge_length = sum(weights .* geom.J_1D)
physical_length = norm([geom.x_phys[end], geom.y_phys[end]] -
                       [geom.x_phys[1], geom.y_phys[1]])
@assert edge_length ≈ physical_length
```
"""
function face_geometry(face::Face, nodes::AbstractVector, corners::AbstractMatrix,
                      jac_matrix::AbstractArray)
    nnodes = length(nodes)

    # Get reference coordinates for this face
    ξ_vals, η_vals = reference_coordinates(face.local_face_id, nodes)

    # Apply flip if needed
    if face.flip
        ξ_vals = reverse(ξ_vals)
        η_vals = reverse(η_vals)
    end

    # Allocate output
    x_phys = zeros(nnodes)
    y_phys = zeros(nnodes)
    t_x = zeros(nnodes)
    t_y = zeros(nnodes)
    n_x = zeros(nnodes)
    n_y = zeros(nnodes)
    J_1D = zeros(nnodes)

    # Get indices for accessing Jacobian matrix
    i_idx, j_idx = canonical_face_nodes(face.local_face_id, nnodes)

    # For each node along the face
    for (idx, (ξ, η)) in enumerate(zip(ξ_vals, η_vals))
        # Physical coordinates (could also get from node_coordinates array)
        x_phys[idx], y_phys[idx] = straight_side_quad_map(ξ, η, corners)

        # Get Jacobian matrix at this point
        # jac_matrix[spatial_dim, ref_dim, i, j]
        # J = [∂x/∂ξ  ∂x/∂η]
        #     [∂y/∂ξ  ∂y/∂η]

        # Determine which node index this corresponds to
        # NOTE: When face.flip=true, we've already reversed ξ_vals/η_vals,
        # but the Jacobian matrix node indices go in canonical order.
        # So we need to reverse idx when flipped.
        node_idx = face.flip ? (nnodes - idx + 1) : idx

        if face.local_face_id == 1  # Bottom: j=1, i varies
            i_node = node_idx
            j_node = 1
        elseif face.local_face_id == 2  # Right: i=nnodes, j varies
            i_node = nnodes
            j_node = node_idx
        elseif face.local_face_id == 3  # Top: j=nnodes, i varies
            i_node = node_idx
            j_node = nnodes
        elseif face.local_face_id == 4  # Left: i=1, j varies
            i_node = 1
            j_node = node_idx
        end

        # Extract Jacobian components
        ∂x_∂ξ = jac_matrix[1, 1, i_node, j_node]
        ∂x_∂η = jac_matrix[1, 2, i_node, j_node]
        ∂y_∂ξ = jac_matrix[2, 1, i_node, j_node]
        ∂y_∂η = jac_matrix[2, 2, i_node, j_node]

        # Compute tangent vector along the face
        # Face 1,3 (bottom/top): s = ξ, so t = ∂x/∂ξ
        # Face 2,4 (right/left): s = η, so t = ∂x/∂η
        if face.local_face_id == 1 || face.local_face_id == 3
            t_x[idx] = ∂x_∂ξ
            t_y[idx] = ∂y_∂ξ
        else  # face 2 or 4
            t_x[idx] = ∂x_∂η
            t_y[idx] = ∂y_∂η
        end

        # Line Jacobian: |t|
        J_1D[idx] = sqrt(t_x[idx]^2 + t_y[idx]^2)

        # Outward normal (2D): rotate tangent by 90°
        # For bottom/top faces (η constant): n = (-∂y/∂ξ, ∂x/∂ξ) or (∂y/∂ξ, -∂x/∂ξ)
        # For left/right faces (ξ constant): n = (∂y/∂η, -∂x/∂η) or (-∂y/∂η, ∂x/∂η)
        #
        # General rule: n = (t_y, -t_x) for outward pointing
        # But need to check sign based on face orientation

        # Compute 2D Jacobian determinant for sign
        Jdet = ∂x_∂ξ * ∂y_∂η - ∂x_∂η * ∂y_∂ξ
        s = sign(Jdet)

        # Outward normals (following Geometry.jl convention)
        if face.local_face_id == 1  # Bottom
            n_x[idx] = -s * (-∂y_∂ξ)
            n_y[idx] = -s * ∂x_∂ξ
        elseif face.local_face_id == 2  # Right
            n_x[idx] = s * ∂y_∂η
            n_y[idx] = s * (-∂x_∂η)
        elseif face.local_face_id == 3  # Top
            n_x[idx] = s * (-∂y_∂ξ)
            n_y[idx] = s * ∂x_∂ξ
        elseif face.local_face_id == 4  # Left
            n_x[idx] = -s * ∂y_∂η
            n_y[idx] = -s * (-∂x_∂η)
        end
    end

    return (
        x_phys = x_phys,
        y_phys = y_phys,
        t_x = t_x,
        t_y = t_y,
        n_x = n_x,
        n_y = n_y,
        J_1D = J_1D
    )
end

"""
    build_face_map(mesh::UnstructuredMesh2D) -> FaceMap

Build complete face information table from mesh connectivity.

# Arguments
- `mesh::UnstructuredMesh2D`: Trixi.jl mesh

# Returns
- `FaceMap`: Complete face table with interior + boundary faces

# Algorithm
1. Loop through mesh.neighbour_information
2. For each entry, create FaceInfo:
   - If neighbor_element == 0: boundary face
   - Otherwise: interior face with neighbor info
3. Determine flip bit by comparing node orderings
4. Group boundary faces by name

# Notes
This builds the "single source of truth" for all face operations.
Call once at mesh setup, then query for all boundary/flux operations.

# Example
```julia
face_map = build_face_map(mesh)

# Get all fault boundary faces
fault_faces = [face_map.faces[i] for i in face_map.boundary_faces[:fault]]

# Get geometry for first fault face
geom = face_geometry(fault_faces[1].face, nodes, corners, jac_matrix)
```
"""
function build_face_map(mesh::UnstructuredMesh2D)
    faces = FaceInfo[]
    boundary_faces = Dict{Symbol, Vector{Int}}()

    # mesh.neighbour_information structure:
    # [1, i] = node1_idx (corner)
    # [2, i] = node2_idx (corner)
    # [3, i] = element1
    # [4, i] = element2 (0 for boundary)
    # [5, i] = surface1 (local face ID for element1)
    # [6, i] = surface2 (local face ID for element2, 0 for boundary)

    for i in 1:size(mesh.neighbour_information, 2)
        el1 = mesh.neighbour_information[3, i]
        el2 = mesh.neighbour_information[4, i]
        surf1_raw = mesh.neighbour_information[5, i]
        surf2_raw = mesh.neighbour_information[6, i]

        # Extract canonical face ID (1-4) and initial flip from sign
        # Negative surface ID means reversed orientation in Trixi.jl
        surf1 = abs(surf1_raw)
        flip1 = surf1_raw < 0

        # Boundary face?
        if el2 == 0 || surf2_raw == 0
            # Get boundary name from mesh
            # mesh.boundary_names is indexed as: (element-1)*4 + surface_id
            # where surface_id is the absolute value (canonical 1-4)
            linear_idx = (el1 - 1) * 4 + surf1
            boundary_name = mesh.boundary_names[linear_idx]

            # Create face with initial flip from mesh orientation
            face = Face(el1, surf1, flip1)
            face_info = FaceInfo(face, boundary_name)

            push!(faces, face_info)

            # Add to boundary map
            if !haskey(boundary_faces, boundary_name)
                boundary_faces[boundary_name] = Int[]
            end
            push!(boundary_faces[boundary_name], length(faces))

        else
            # Interior face
            surf2 = abs(surf2_raw)
            flip2 = surf2_raw < 0

            # Initial flip from mesh orientation
            # This will be refined in determine_face_flips! based on DOF connectivity
            face = Face(el1, surf1, flip1)
            face_info = FaceInfo(face, el2, surf2, flip2)

            push!(faces, face_info)
        end
    end

    return FaceMap(faces, boundary_faces)
end

"""
    determine_face_flips!(face_map::FaceMap, dof_id::Array{Int,3})

Determine flip bits for all faces based on DOF connectivity.

# Arguments
- `face_map::FaceMap`: Face map (modified in-place)
- `dof_id::Array{Int,3}`: Connectivity matrix [nnodes, nnodes, nelements]

# Notes
For each interior face, compare node orderings:
- If face1_nodes == face2_nodes: no flip
- If face1_nodes == reverse(face2_nodes): flip = true

Modifies face_map.faces in place to set correct flip bits.

# Example
```julia
face_map = build_face_map(mesh)
determine_face_flips!(face_map, dof_id)
```
"""
function determine_face_flips!(face_map::FaceMap, dof_id::Array{Int,3})
    nnodes = size(dof_id, 1)

    for (idx, face_info) in enumerate(face_map.faces)
        # Skip boundary faces
        if face_info.type == Boundary
            continue
        end

        # Get node IDs for this face
        face1 = face_info.face
        el1 = face1.element
        surf1 = face1.local_face_id

        # Get neighbor
        el2 = face_info.neighbor_element
        surf2 = face_info.neighbor_local_face_id

        # Extract node IDs (without flip first)
        nodes1 = extract_face_values(dof_id[:, :, el1], Face(el1, surf1, false))
        nodes2 = extract_face_values(dof_id[:, :, el2], Face(el2, surf2, false))

        # Determine if reversed
        if nodes1 == nodes2
            # Same ordering, no flip
            flip = false
        elseif nodes1 == reverse(nodes2)
            # Reversed ordering, flip = true
            flip = true
        else
            @warn "Face nodes don't match!" el1 surf1 el2 surf2 nodes1 nodes2
            flip = false
        end

        # Update face with correct flip
        # (Julia structs are immutable, need to reconstruct)
        new_face = Face(el1, surf1, flip)
        new_face_info = FaceInfo(new_face, el2, surf2, false)  # neighbor flip TBD

        face_map.faces[idx] = new_face_info
    end

    return face_map
end

# Forward declaration - this function should be provided by including module
# Either import from Geometry.jl or define before calling face_geometry()
if !@isdefined(straight_side_quad_map)
    """
    Placeholder for straight_side_quad_map - should be defined by including module.
    See Geometry.jl for implementation.
    """
    function straight_side_quad_map end
end
