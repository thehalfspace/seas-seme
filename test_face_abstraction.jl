"""
Test script to validate the new Face abstraction against the old implementation.

This script:
1. Loads both old and new boundary extraction code
2. Builds a test mesh
3. Compares results from both implementations
4. Runs geometric invariant checks
5. Reports any discrepancies

Run with geometric checks:
    SEAS_SEME_GEOMETRIC_CHECKS=1 julia test_face_abstraction.jl
"""

using Trixi
using FastGaussQuadrature

# Add source to path
push!(LOAD_PATH, joinpath(@__DIR__, "src"))

# Load mesh modules
include("src/Mesh/Geometry.jl")
include("src/Mesh/Connectivity.jl")
include("src/Mesh/Boundaries.jl")  # Old implementation
include("src/Mesh/BoundariesNew.jl")  # New implementation

# Make Geometry functions available
using .SEAS_SEME.Mesh: straight_side_quad_map

println("="^80)
println("Testing Face Abstraction vs Old Implementation")
println("="^80)

# Load test mesh
mesh_file = "assets/mesh_2d_homogeneous.mesh"
if !isfile(mesh_file)
    error("Mesh file not found: $mesh_file")
end

println("\n📂 Loading mesh: $mesh_file")
mesh = UnstructuredMesh2D(mesh_file)

println("   Elements: $(mesh.n_elements)")
println("   Polynomial degree: $(mesh.polydeg)")
println("   Boundaries: $(unique(mesh.boundary_names))")

# Setup
polydeg = mesh.polydeg
nnodes = polydeg + 1
nel = mesh.n_elements

# Get GLL nodes and weights
nodes, weights = gausslobatto(nnodes)

# Build geometry
println("\n🔧 Building mesh geometry...")

# Node coordinates
node_coordinates = zeros(2, nnodes, nnodes, nel)
for el in 1:nel
    corners = mesh.corners[:, mesh.element_node_ids[:, el]]'
    calc_node_coordinates!(node_coordinates, el, nodes, corners)
end

# Jacobian matrices
jacobian_matrix = zeros(2, 2, nnodes, nnodes, nel)
for el in 1:nel
    corners = mesh.corners[:, mesh.element_node_ids[:, el]]'
    calc_metric_terms!(jacobian_matrix, el, nodes, corners)
end

# Connectivity
println("🔗 Building connectivity...")
dof_id = connectivity_matrix(mesh)
n_dofs = maximum(dof_id)
println("   Total DOFs: $n_dofs")

# Test parameters
impedance_absorbing = 1.0  # Would be ρ*vs in real simulation
impedance_fault = 1.0

# Test boundaries
boundaries_to_test = [:fault, :creep, :absorbing, :free_surface]

println("\n" * "="^80)
println("Running Comparison Tests")
println("="^80)

all_match = true

for boundary_name in boundaries_to_test
    println("\n🧪 Testing boundary: $boundary_name")

    impedance = (boundary_name == :absorbing) ? impedance_absorbing : impedance_fault

    try
        old_result, new_result, match = compare_boundary_implementations(
            mesh, node_coordinates, jacobian_matrix, weights,
            impedance, dof_id, boundary_name
        )

        old_ids, old_x, old_y, old_mat = old_result
        new_ids, new_x, new_y, new_mat = new_result

        if match
            println("   ✅ PASS - Results identical!")
            println("      Nodes: $(length(old_ids))")
            println("      Max impedance: $(maximum(old_mat))")
        else
            println("   ❌ FAIL - Results differ!")
            all_match = false

            # Show detailed differences
            if old_ids != new_ids
                println("      Node ID mismatch:")
                println("        Old: $(length(old_ids)) nodes")
                println("        New: $(length(new_ids)) nodes")
            end

            if !isapprox(old_x, new_x)
                println("      X coordinate mismatch:")
                println("        Max error: $(maximum(abs.(old_x .- new_x)))")
            end

            if !isapprox(old_y, new_y)
                println("      Y coordinate mismatch:")
                println("        Max error: $(maximum(abs.(old_y .- new_y)))")
            end

            if !isapprox(old_mat, new_mat)
                println("      Impedance matrix mismatch:")
                println("        Max error: $(maximum(abs.(old_mat .- new_mat)))")
                println("        Relative error: $(maximum(abs.(old_mat .- new_mat) ./ (old_mat .+ 1e-10)))")
            end
        end

    catch e
        println("   ❌ ERROR: $e")
        all_match = false
        showerror(stdout, e, catch_backtrace())
        println()
    end
end

# Test geometric checks (if enabled)
if get(ENV, "SEAS_SEME_GEOMETRIC_CHECKS", "0") == "1"
    println("\n" * "="^80)
    println("Running Geometric Invariant Checks")
    println("="^80)

    println("\n🔍 Checking Jacobian determinants...")
    try
        check_jacobian_determinant_all(jacobian_matrix)
        println("   ✅ All Jacobian determinants positive")
    catch e
        println("   ❌ Jacobian check failed: $e")
    end

    println("\n🔍 Checking edge lengths...")
    # Build face map
    face_map = build_face_map(mesh)
    determine_face_flips!(face_map, dof_id)

    edge_check_passed = true
    for (boundary_name, face_indices) in face_map.boundary_faces
        for idx in face_indices
            face_info = face_map.faces[idx]
            face = face_info.face
            el = face.element

            corners = mesh.corners[:, mesh.element_node_ids[:, el]]'
            geom = face_geometry(face, nodes, corners, jacobian_matrix[:, :, :, :, el])

            try
                check_edge_length(geom, weights)
            catch e
                println("   ❌ Edge length check failed: boundary=$boundary_name, el=$el, face=$(face.local_face_id)")
                edge_check_passed = false
            end
        end
    end

    if edge_check_passed
        println("   ✅ All edge lengths consistent")
    end
end

# Summary
println("\n" * "="^80)
println("Test Summary")
println("="^80)

if all_match
    println("✅ ALL TESTS PASSED")
    println("\nThe new Face abstraction produces identical results to the old implementation.")
    println("Safe to migrate to the new unbreakable version!")
else
    println("❌ SOME TESTS FAILED")
    println("\nInvestigate discrepancies before migrating.")
end

println("\n" * "="^80)
