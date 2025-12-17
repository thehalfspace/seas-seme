"""
Quick integration test: build_mesh with Face abstraction integrated.
"""

using Trixi

println("="^80)
println("Integration Test: build_mesh with Face Abstraction")
println("="^80)

# Load mesh modules in correct order
include("src/Config/Parameters.jl")
include("src/Mesh/Connectivity.jl")
include("src/Mesh/Geometry.jl")
include("src/Mesh/Faces.jl")
include("src/Mesh/GeometricChecks.jl")
include("src/Mesh/Boundaries.jl")
include("src/Mesh/Mesh.jl")

# Create test config
struct TestPhysicsConfig
    density::Float64
    shear_velocity::Float64
end

struct TestMeshConfig
    file::String
    polynomial_degree::Int
end

physics = TestPhysicsConfig(2670.0, 3464.0)
mesh_config = TestMeshConfig("data/mesh/unstructured/unstructured.mesh", 4)

println("\n🧪 Testing build_mesh...")
try
    mesh = build_mesh(mesh_config, physics)

    println("\n✅ Mesh built successfully!")
    println("   Elements: $(mesh.n_elements)")
    println("   DOFs: $(mesh.ndof)")
    println("   Polynomial degree: $(mesh.polynomial_degree)")
    println("   Face map faces: $(length(mesh.face_map.faces))")
    println("   Boundary faces: $(keys(mesh.face_map.boundary_faces))")

    # Verify face_map is populated
    for (boundary_name, face_indices) in mesh.face_map.boundary_faces
        println("   $boundary_name: $(length(face_indices)) faces")
    end

    println("\n✅ INTEGRATION TEST PASSED")

catch e
    println("\n❌ INTEGRATION TEST FAILED")
    println("Error: $e")
    showerror(stdout, e, catch_backtrace())
    println()
end

println("\n" * "="^80)
