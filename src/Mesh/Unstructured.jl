"""
Unstructured mesh generation using HOHQMesh.jl

This module provides functionality to generate unstructured quadrilateral meshes
for SEAS-SEME earthquake cycle simulations using the HOHQMesh.jl package.

# Mesh Geometry

Standard fault mesh configuration:

```
           Free Surface
    ++++++++++++++++++++++++++++++++++(Lx, Ly)
    |                                :
    |                                :
    |                                :  Fault Boundary
    |                                :
    |                                :
    |                                (Lx, fault_y)
    |                                }
    |                                }   Creep Boundary
    |                                }
    |                                }
(0,0)-------------------------------(Lx, 0)
       Absorbing boundaries: left and bottom
```

# Boundary Definitions
- **Fault**: Right edge from fault_y to Ly (where earthquakes nucleate)
- **Creep**: Right edge from 0 to fault_y (constant velocity plate loading)
- **Absorbing**: Left edge and bottom edge (prevent wave reflections)
- **Free Surface**: Top edge (stress-free boundary)

# References
- HOHQMesh documentation: https://trixi-framework.github.io/HOHQMesh/
- Trixi.jl unstructured meshes: https://trixi-framework.github.io/Trixi.jl/stable/meshes/unstructured_quad_mesh/
"""

using HOHQMesh

"""
    MeshParameters

Configuration parameters for unstructured mesh generation.

# Fields
- `Lx::Float64`: Domain width (meters)
- `Ly::Float64`: Domain height (meters)
- `fault_start_y::Float64`: Y-coordinate where fault begins (meters)
- `background_grid_size::Float64`: Default element size (meters)
- `fault_refinement_h::Float64`: Target element size along fault (meters)
- `fault_refinement_width::Float64`: Width of refinement zone around fault (meters)
- `polynomial_order::Int`: Polynomial degree for spectral elements (default: 4)
- `smoothing_iterations::Int`: Number of mesh smoothing iterations (default: 50)
- `dip_angle::Float64`: Fault dip angle in degrees (90 = vertical, default: 90.0)
"""
struct MeshParameters
    Lx::Float64
    Ly::Float64
    fault_start_y::Float64
    background_grid_size::Float64
    fault_refinement_h::Float64
    fault_refinement_width::Float64
    polynomial_order::Int
    smoothing_iterations::Int
    dip_angle::Float64
end

"""
    MeshParameters(;kwargs...)

Construct mesh parameters with keyword arguments.

# Keyword Arguments
- `Lx=40.0e3`: Domain width (meters)
- `Ly=40.0e3`: Domain height (meters)
- `fault_start_y=20.0e3`: Y-coordinate where fault begins (meters)
- `background_grid_size=2.0e3`: Default element size (meters)
- `fault_refinement_h=0.2e3`: Target element size along fault (meters)
- `fault_refinement_width=2.0e3`: Width of refinement zone (meters)
- `polynomial_order=4`: Polynomial degree
- `smoothing_iterations=50`: Smoothing iterations
- `dip_angle=90.0`: Fault dip angle in degrees (90 = vertical)

# Example
```julia
params = MeshParameters(
    Lx=40.0e3,
    Ly=40.0e3,
    fault_start_y=20.0e3,
    background_grid_size=2.0e3,
    fault_refinement_h=0.2e3,
    dip_angle=60.0
)
```
"""
function MeshParameters(;
    Lx::Float64 = 40.0e3,
    Ly::Float64 = 40.0e3,
    fault_start_y::Float64 = 20.0e3,
    background_grid_size::Float64 = 2.0e3,
    fault_refinement_h::Float64 = 0.2e3,
    fault_refinement_width::Float64 = 2.0e3,
    polynomial_order::Int = 4,
    smoothing_iterations::Int = 50,
    dip_angle::Float64 = 90.0
)
    MeshParameters(
        Lx, Ly, fault_start_y,
        background_grid_size,
        fault_refinement_h,
        fault_refinement_width,
        polynomial_order,
        smoothing_iterations,
        dip_angle
    )
end

"""
    create_control_file(params::MeshParameters, control_file_path::String)

Generate HOHQMesh control file for unstructured mesh.

# Arguments
- `params::MeshParameters`: Mesh generation parameters
- `control_file_path::String`: Path where control file will be written

# File Format
The control file uses HOHQMesh's specification format with sections:
- CONTROL_INPUT: Run parameters and mesh settings
- BACKGROUND_GRID: Default element sizing
- SPRING_SMOOTHER: Mesh smoothing configuration
- REFINEMENT_REGIONS: Fault refinement zone
- MODEL: Boundary definitions (fault, creep, absorbing, free surface)

# Example
```julia
params = MeshParameters()
create_control_file(params, "mesh.control")
```
"""
function create_control_file(params::MeshParameters, control_file_path::String)
    δ = params.dip_angle
    Lx = params.Lx
    Ly = params.Ly
    fy = params.fault_start_y

    # Compute boundary vertex coordinates for dipping fault geometry
    # For δ=90°: x_b = Lx, x_t = Lx (rectangular domain, recovers current geometry)
    # For δ<90°: domain becomes a pentagon with dipping right edge
    if δ == 90.0
        x_b = Lx
        x_t = Lx
    else
        δ_rad = deg2rad(δ)
        x_b = Lx - Ly / tan(δ_rad)
        x_t = Lx - (Ly - fy) / tan(δ_rad)
        x_b > 0 || error("Domain too narrow for dip angle $(δ)°: need Lx > Ly/tan(δ) = $(Ly/tan(δ_rad)/1e3) km, got Lx = $(Lx/1e3) km")
    end

    open(control_file_path, "w") do io
        # Write control file in the exact format that works
        println(io, "\t\t\\begin{CONTROL_INPUT}")
        println(io, "\t\t\t\\begin{RUN_PARAMETERS}")
        println(io, "\t\t\t\tmesh file name   = unstructured.mesh")
        println(io, "\t\t\t\tplot file name   = unstructured.tec")
        println(io, "\t\t\t\tstats file name  = none")
        println(io, "\t\t\t\tmesh file format = ISM-v2")
        println(io, "\t\t\t\tpolynomial order = $(params.polynomial_order)")
        println(io, "\t\t\t\tplot file format = skeleton")
        println(io, "\t\t\t\\end{RUN_PARAMETERS}")
        println(io, "")
        println(io, "\t\t\t\\begin{BACKGROUND_GRID}")
        println(io, "\t\t\t\tbackground grid size = [$(params.background_grid_size),$(params.background_grid_size),0.0]")
        println(io, "\t\t\t\\end{BACKGROUND_GRID}")
        println(io, "")
        println(io, "\t\t\t\\begin{SPRING_SMOOTHER}")
        println(io, "\t\t\t\tsmoothing            = ON")
        println(io, "\t\t\t\tsmoothing type       = LinearAndCrossBarSpring")
        println(io, "\t\t\t\tnumber of iterations = $(params.smoothing_iterations)")
        println(io, "\t\t\t\\end{SPRING_SMOOTHER}")
        println(io, "\t\t\t\\begin{REFINEMENT_REGIONS}")
        println(io, "\t\t\t\\begin{REFINEMENT_LINE}")
        println(io, "\t\t\t\ttype = smooth")
        println(io, "\t\t\t\tx0   = [$(x_t),$(fy),0.0]")
        println(io, "\t\t\t\tx1   = [$(Lx),$(Ly),0.0]")
        println(io, "\t\t\t\th    = $(params.fault_refinement_h)")
        println(io, "\t\t\t\tw    = $(params.fault_refinement_width)")
        println(io, "\t\t\t\\end{REFINEMENT_LINE}")
        println(io, "\t\t\t\\end{REFINEMENT_REGIONS}")
        println(io, "\t\t\\end{CONTROL_INPUT}")
        println(io, "\t\t")
        println(io, "\t\t\\begin{MODEL}")
        println(io, "\t\t\t\\begin{OUTER_BOUNDARY}")

        # Segment 1: bottom absorbing — (0,0) → (x_b, 0)
        println(io, "\t\t\t\\begin{END_POINTS_LINE}")
        println(io, "\t\t\t\tname = absorbing ")
        println(io, "\t\t\t\txStart = [0.0,0.0,0.0]")
        println(io, "\t\t\t\txEnd   = [$(x_b),0.0,0.0]")
        println(io, "\t\t\t\\end{END_POINTS_LINE}")

        # Segment 2: creep — (x_b, 0) → (x_t, fy)
        println(io, "\t\t\t\\begin{END_POINTS_LINE}")
        println(io, "\t\t\t\tname = creep")
        println(io, "\t\t\t\txStart = [$(x_b),0.0,0.0]")
        println(io, "\t\t\t\txEnd   = [$(x_t),$(fy),0.0]")
        println(io, "\t\t\t\\end{END_POINTS_LINE}")

        # Segment 3: fault — (x_t, fy) → (Lx, Ly)
        println(io, "\t\t\t\\begin{END_POINTS_LINE}")
        println(io, "\t\t\t\tname = fault")
        println(io, "\t\t\t\txStart = [$(x_t),$(fy),0.0]")
        println(io, "\t\t\t\txEnd   = [$(Lx),$(Ly),0.0]")
        println(io, "\t\t\t\\end{END_POINTS_LINE}")

        # Segment 4: free surface — (Lx, Ly) → (0, Ly)
        println(io, "\t\t\t\\begin{END_POINTS_LINE}")
        println(io, "\t\t\t\tname = free_surface")
        println(io, "\t\t\t\txStart = [$(Lx),$(Ly),0.0]")
        println(io, "\t\t\t\txEnd   = [0.0,$(Ly),0.0]")
        println(io, "\t\t\t\\end{END_POINTS_LINE}")

        # Segment 5: left absorbing — (0, Ly) → (0, 0)
        println(io, "\t\t\t\\begin{END_POINTS_LINE}")
        println(io, "\t\t\t\tname = absorbing ")
        println(io, "\t\t\t\txStart = [0.0,$(Ly),0.0]")
        println(io, "\t\t\t\txEnd   = [0.0,0.0,0.0]")
        println(io, "\t\t\t\\end{END_POINTS_LINE}")

        println(io, "\t\t\t\\end{OUTER_BOUNDARY}")
        println(io, "\t\t\\end{MODEL}")
        println(io, "\t\t\\end{FILE}")
    end
end

"""
    generate_unstructured_mesh(params::MeshParameters, output_dir::String;
                               mesh_name::String="unstructured")

Generate unstructured quadrilateral mesh using HOHQMesh.

# Arguments
- `params::MeshParameters`: Mesh generation parameters
- `output_dir::String`: Directory where mesh files will be saved

# Keyword Arguments
- `mesh_name::String="unstructured"`: Base name for mesh files

# Output Files
The function generates several files in `output_dir`:
- `{mesh_name}.control`: HOHQMesh control file
- `{mesh_name}.mesh`: Binary mesh file (ISM-v2 format)
- `{mesh_name}.tec`: Tecplot skeleton file for visualization

# Example
```julia
params = MeshParameters(
    Lx=40.0e3,
    Ly=40.0e3,
    fault_start_y=20.0e3
)
generate_unstructured_mesh(params, "data/mesh/my_simulation")
```

# Notes
- Creates output directory if it doesn't exist
- Overwrites existing mesh files in output_dir
- Use ParaView to visualize the .tec file
"""
function generate_unstructured_mesh(
    params::MeshParameters,
    output_dir::String;
    mesh_name::String = "unstructured"
)
    # Create output directory
    if !isdir(output_dir)
        mkpath(output_dir)
        @info "Created mesh directory: $output_dir"
    else
        @info "Using existing mesh directory: $output_dir"
    end

    # Generate control file
    control_file = joinpath(output_dir, "$(mesh_name).control")
    @info "Writing control file: $control_file"
    create_control_file(params, control_file)

    # Generate mesh
    @info "Generating mesh with HOHQMesh..."
    @info "  Domain: $(params.Lx/1e3) × $(params.Ly/1e3) km"
    @info "  Fault location: y = $(params.fault_start_y/1e3) km"
    @info "  Dip angle: $(params.dip_angle)°"
    @info "  Background element size: $(params.background_grid_size/1e3) km"
    @info "  Fault refinement: $(params.fault_refinement_h/1e3) km"
    @info "  Polynomial order: $(params.polynomial_order)"

    generate_mesh(control_file, output_directory=output_dir)

    mesh_file = joinpath(output_dir, "$(mesh_name).mesh")
    if isfile(mesh_file)
        @info "✓ Mesh generation complete: $mesh_file"
        @info "  Visualize with ParaView: $(joinpath(output_dir, "$(mesh_name).tec"))"
        return mesh_file
    else
        error("Mesh generation failed - mesh file not created")
    end
end

"""
    generate_default_mesh(output_dir::String)

Generate default unstructured mesh with standard parameters.

# Arguments
- `output_dir::String`: Directory where mesh files will be saved

# Default Parameters
- Domain: 40 × 40 km
- Fault start: 20 km
- Background element size: 2 km
- Fault refinement: 0.2 km
- Polynomial order: 4

# Example
```julia
mesh_file = generate_default_mesh("data/mesh/strike_slip_benchmark")
```
"""
function generate_default_mesh(output_dir::String)
    params = MeshParameters()
    return generate_unstructured_mesh(params, output_dir)
end
