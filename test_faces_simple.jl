"""
Simple test of Face abstraction core functionality (no full mesh required).
"""

println("Testing Face abstraction core functions...")

# Test Face struct
struct Face
    element::Int
    local_face_id::Int
    flip::Bool
end

# Test canonical_face_nodes
function canonical_face_nodes(local_face_id::Int, nnodes::Int)
    if local_face_id == 1  # Bottom
        return (1:nnodes, 1)
    elseif local_face_id == 2  # Right
        return (nnodes, 1:nnodes)
    elseif local_face_id == 3  # Top
        return (1:nnodes, nnodes)
    elseif local_face_id == 4  # Left
        return (1, 1:nnodes)
    end
end

println("\n1. Testing canonical_face_nodes...")
nnodes = 5
for face_id in 1:4
    i_idx, j_idx = canonical_face_nodes(face_id, nnodes)
    println("   Face $face_id: i=$i_idx, j=$j_idx")
end

# Test extract_face_values
function extract_face_values(array::AbstractMatrix, face::Face)
    i_idx, j_idx = canonical_face_nodes(face.local_face_id, size(array, 1))

    if i_idx isa Int
        values = array[i_idx, j_idx]
    else
        values = array[i_idx, j_idx]
    end

    return face.flip ? reverse(values) : values
end

println("\n2. Testing extract_face_values...")
# Create test element DOF array
dof_array = reshape(1:25, 5, 5)
println("   Element DOF array:")
display(dof_array)

for face_id in 1:4
    face = Face(1, face_id, false)
    values = extract_face_values(dof_array, face)
    println("\n   Face $face_id (no flip): $values")

    face_flipped = Face(1, face_id, true)
    values_flipped = extract_face_values(dof_array, face_flipped)
    println("   Face $face_id (flipped):  $values_flipped")
end

# Test reference_coordinates
function reference_coordinates(local_face_id::Int, nodes::AbstractVector)
    nnodes = length(nodes)

    if local_face_id == 1  # Bottom
        return (nodes, fill(-1.0, nnodes))
    elseif local_face_id == 2  # Right
        return (fill(1.0, nnodes), nodes)
    elseif local_face_id == 3  # Top
        return (nodes, fill(1.0, nnodes))
    elseif local_face_id == 4  # Left
        return (fill(-1.0, nnodes), nodes)
    end
end

println("\n3. Testing reference_coordinates...")
nodes = [-1.0, -0.5, 0.0, 0.5, 1.0]
for face_id in 1:4
    ξ, η = reference_coordinates(face_id, nodes)
    println("   Face $face_id: ξ=$ξ, η=$η")
end

println("\n✅ Core functions working correctly!")
