"""
Validate face_geometry against current Boundaries.jl implementation.

This test compares the Jacobian extraction from:
1. Old implementation: Manual slicing (jac[1,1,:,1], jac[1,2,end,:], etc.)
2. New implementation: face_geometry() abstraction

Should produce identical J_1D values.
"""

using Trixi
using LinearAlgebra

println("="^80)
println("Validating face_geometry vs Boundaries.jl")
println("="^80)

# Load implementations
include("src/Mesh/Geometry.jl")
include("src/Mesh/Connectivity.jl")
include("src/Mesh/Boundaries.jl")  # Old implementation
include("src/Mesh/Faces.jl")       # New implementation

# Load test mesh (unstructured)
mesh_file = "data/mesh/unstructured/unstructured.mesh"
if !isfile(mesh_file)
    error("Mesh file not found: $mesh_file")
end

println("\n📂 Loading mesh: $mesh_file")
mesh = UnstructuredMesh2D(mesh_file)

println("   Elements: $(mesh.n_elements)")
println("   Polynomial degree: $(mesh.polydeg)")
println("   Boundaries: $(unique(mesh.boundary_names))")

# Setup geometry
polydeg = mesh.polydeg
nnodes = polydeg + 1
nel = mesh.n_elements

# Get GLL nodes via Trixi's LobattoLegendreBasis
basis = LobattoLegendreBasis(polydeg)
nodes = basis.nodes
weights = basis.weights

println("\n🔧 Building geometry...")

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

# Build connectivity
println("🔗 Building connectivity...")
dof_id = connectivity_matrix(mesh)

println("\n" * "="^80)
println("Comparing J_1D Computation Methods")
println("="^80)

# Test all four boundary types
boundaries_to_test = [:fault, :creep, :absorbing, :free_surface]

global all_match = true
global max_error_overall = 0.0

for boundary_name in boundaries_to_test
    println("\n🧪 Testing boundary: $boundary_name")

    # Get boundary elements from old implementation
    boundary_el_id = if boundary_name == :absorbing
        findall(mesh.boundary_names .== :absorbing)
    elseif boundary_name == :free_surface
        findall(mesh.boundary_names .== :free_surface)
    elseif boundary_name == :fault
        findall(mesh.boundary_names .== :fault)
    elseif boundary_name == :creep
        findall(mesh.boundary_names .== :creep)
    end

    if isempty(boundary_el_id)
        println("   ⚠️  Boundary not found in mesh, skipping")
        continue
    end

    println("   Found $(length(boundary_el_id)) boundary faces")

    # Compare J_1D for each boundary face
    max_error = 0.0
    match_count = 0
    mismatch_count = 0

    for id in boundary_el_id
        surface = id[1]
        el = id[2]

        corners = mesh.corners[:, mesh.element_node_ids[:, el]]'
        jac = jacobian_matrix[:, :, :, :, el]

        # OLD METHOD: Manual Jacobian slicing (from Boundaries.jl:193-207)
        jac1D_old = if surface == 1  # Bottom
            sqrt.(jac[1, 1, :, 1] .^ 2 .+ jac[2, 1, :, 1] .^ 2)
        elseif surface == 2  # Right
            sqrt.(jac[1, 2, end, :] .^ 2 .+ jac[2, 2, end, :] .^ 2)
        elseif surface == 3  # Top
            sqrt.(jac[1, 1, :, end] .^ 2 .+ jac[2, 1, :, end] .^ 2)
        elseif surface == 4  # Left
            sqrt.(jac[1, 2, 1, :] .^ 2 .+ jac[2, 2, 1, :] .^ 2)
        else
            error("Invalid surface: $surface")
        end

        # NEW METHOD: face_geometry abstraction
        face = Face(el, surface, false)  # No flip for boundary faces initially
        geom = face_geometry(face, nodes, corners, jac)
        jac1D_new = geom.J_1D

        # Compare
        error = maximum(abs.(jac1D_old .- jac1D_new))
        max_error = max(max_error, error)
        global max_error_overall = max(max_error_overall, error)

        if error < 1e-10
            match_count += 1
        else
            mismatch_count += 1
            if mismatch_count <= 3  # Show first 3 mismatches
                println("      ❌ Mismatch on element $el, surface $surface:")
                println("         Old J_1D: $(jac1D_old)")
                println("         New J_1D: $(jac1D_new)")
                println("         Error: $error")
            end
        end
    end

    if mismatch_count == 0
        println("   ✅ PASS - All $(match_count) faces match exactly")
        println("      Max error: $(max_error)")
    else
        println("   ❌ FAIL - $(mismatch_count) mismatches out of $(match_count + mismatch_count)")
        println("      Max error: $(max_error)")
        global all_match = false
    end
end

println("\n" * "="^80)
println("Validation Summary")
println("="^80)

if all_match
    println("✅ VALIDATION PASSED")
    println("\nThe new face_geometry() produces identical J_1D values to the old")
    println("manual slicing method (max error: $(max_error_overall)).")
    println("\nSafe to replace Boundaries.jl with the new implementation!")
else
    println("❌ VALIDATION FAILED")
    println("\nDiscrepancies found between old and new implementations.")
    println("Max error: $(max_error_overall)")
    println("\nInvestigate before migrating.")
end

println("\n" * "="^80)
