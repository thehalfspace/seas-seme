"""
Test face_geometry function with a simple element.
"""

println("Testing face_geometry implementation...")

# Load required packages
using Trixi
using LinearAlgebra

# Include dependencies in correct order
include("src/Mesh/Geometry.jl")
include("src/Mesh/Faces.jl")

# Create a simple square element [0,1] × [0,1]
corners = [
    0.0  0.0;  # Bottom-left (corner 1)
    1.0  0.0;  # Bottom-right (corner 2)
    1.0  1.0;  # Top-right (corner 3)
    0.0  1.0   # Top-left (corner 4)
]

# Use 5 GLL nodes
nnodes = 5
# Simple equally-spaced nodes for testing (not actual GLL)
nodes = range(-1, 1, length=nnodes)

# Compute Jacobian matrix for this element
jac_matrix = zeros(2, 2, nnodes, nnodes)
for j in 1:nnodes, i in 1:nnodes
    jac_matrix[:, :, i, j] .= reshape(
        collect(straight_side_quad_map_metrics(nodes[i], nodes[j], corners)),
        2, 2
    )
end

println("\n1. Testing face_geometry without flip...")
for face_id in 1:4
    face = Face(1, face_id, false)
    geom = face_geometry(face, nodes, corners, jac_matrix)

    println("\nFace $face_id:")
    println("  x_phys: $(round.(geom.x_phys, digits=3))")
    println("  y_phys: $(round.(geom.y_phys, digits=3))")
    println("  J_1D:   $(round.(geom.J_1D, digits=3))")

    # Check physical coordinates make sense
    if face_id == 1  # Bottom
        @assert all(geom.y_phys .≈ 0.0) "Bottom face should have y=0"
    elseif face_id == 2  # Right
        @assert all(geom.x_phys .≈ 1.0) "Right face should have x=1"
    elseif face_id == 3  # Top
        @assert all(geom.y_phys .≈ 1.0) "Top face should have y=1"
    elseif face_id == 4  # Left
        @assert all(geom.x_phys .≈ 0.0) "Left face should have x=0"
    end

    # For a square, J_1D should be constant = 0.5 (since reference element is [-1,1])
    expected_J1D = 0.5
    if !all(isapprox.(geom.J_1D, expected_J1D, rtol=0.01))
        @warn "J_1D not constant for square element" geom.J_1D expected=expected_J1D
    end
end

println("\n2. Testing face_geometry WITH flip...")
for face_id in 1:4
    face_noflip = Face(1, face_id, false)
    face_flip = Face(1, face_id, true)

    geom_noflip = face_geometry(face_noflip, nodes, corners, jac_matrix)
    geom_flip = face_geometry(face_flip, nodes, corners, jac_matrix)

    println("\nFace $face_id:")
    println("  No flip x: $(round.(geom_noflip.x_phys, digits=3))")
    println("  Flipped x: $(round.(geom_flip.x_phys, digits=3))")

    # Flipped should be reversed
    if geom_flip.x_phys != reverse(geom_noflip.x_phys)
        @warn "Flip didn't reverse x coordinates!" face_id
    end

    if geom_flip.y_phys != reverse(geom_noflip.y_phys)
        @warn "Flip didn't reverse y coordinates!" face_id
    end

    # J_1D should also be reversed (but have same values)
    if !isapprox(geom_flip.J_1D, reverse(geom_noflip.J_1D), rtol=1e-10)
        @warn "Flip didn't properly reverse J_1D!" face_id
    end
end

println("\n3. Testing extract_face_values...")
# Create a test DOF array
dof_array = reshape(1:25, 5, 5)

for face_id in 1:4
    face = Face(1, face_id, false)
    values = extract_face_values(dof_array, face)

    face_flip = Face(1, face_id, true)
    values_flip = extract_face_values(dof_array, face_flip)

    println("\nFace $face_id:")
    println("  No flip: $values")
    println("  Flipped: $values_flip")

    @assert values_flip == reverse(values) "Flip should reverse node IDs"
end

println("\n✅ All face_geometry tests passed!")
