"""
Mesh Generation Script for SEAS-SEME

This script generates unstructured quadrilateral meshes for earthquake cycle
simulations using HOHQMesh.jl.

# Usage

## From configuration file (recommended):
```bash
julia scripts/generate_mesh.jl config/mesh_default.toml simulation_name
```

## From command line with custom parameters:
```bash
julia scripts/generate_mesh.jl simulation_name --Lx 60e3 --fault-y 30e3 --refine 0.1e3
```

# Command Line Arguments

Positional:
- `config_or_name`: Either a .toml config file path OR simulation name

When using config file:
- `simulation_name`: Name of simulation (creates data/mesh/{simulation_name}/)

When using simulation name directly, optional arguments:
- `--Lx`: Domain width in meters (default: 40e3)
- `--Ly`: Domain height in meters (default: 40e3)
- `--fault-y`: Y-coordinate where fault starts in meters (default: 20e3)
- `--background`: Background element size in meters (default: 2e3)
- `--refine`: Fault refinement element size in meters (default: 0.2e3)
- `--refine-width`: Width of refinement zone in meters (default: 2e3)
- `--order`: Polynomial order (default: 4)
- `--smooth`: Number of smoothing iterations (default: 50)

# Examples

## Using configuration file:
```bash
julia scripts/generate_mesh.jl config/mesh_default.toml strike_slip_benchmark
```

## Using command line (for quick tests):
```bash
julia scripts/generate_mesh.jl test_mesh --refine 0.1e3
```

# Output

Creates directory `data/mesh/{simulation_name}/` containing:
- `unstructured.control`: HOHQMesh control file
- `unstructured.mesh`: Binary mesh file (use in config.toml)
- `unstructured.tec`: Tecplot file (visualize in ParaView)

The mesh file path should be referenced in your simulation config.toml:
```toml
[mesh]
file = "data/mesh/{simulation_name}/unstructured.mesh"
polynomial_degree = 4
```
"""

using ArgParse
using TOML
using HOHQMesh

# Include mesh generation functions directly from source
include(joinpath(@__DIR__, "..", "src", "Mesh", "Unstructured.jl"))

function parse_commandline()
    s = ArgParseSettings(
        description = "Generate unstructured mesh for SEAS-SEME earthquake cycle simulation",
        epilog = """
        Examples:
            julia scripts/generate_mesh.jl config/mesh_default.toml my_simulation
            julia scripts/generate_mesh.jl test_mesh --refine 0.1e3
        """
    )

    @add_arg_table! s begin
        "config_or_name"
            help = "Config file (.toml) or simulation name"
            required = true

        "simulation_name"
            help = "Simulation name (required if first arg is config file)"
            required = false

        "--Lx"
            help = "Domain width in meters (ignored if using config)"
            arg_type = Float64
            default = 40.0e3

        "--Ly"
            help = "Domain height in meters (ignored if using config)"
            arg_type = Float64
            default = 40.0e3

        "--fault-y"
            help = "Y-coordinate where fault starts (meters, ignored if using config)"
            arg_type = Float64
            default = 20.0e3

        "--background"
            help = "Background element size (meters, ignored if using config)"
            arg_type = Float64
            default = 2.0e3

        "--refine"
            help = "Fault refinement element size (meters, ignored if using config)"
            arg_type = Float64
            default = 0.2e3

        "--refine-width"
            help = "Width of fault refinement zone (meters, ignored if using config)"
            arg_type = Float64
            default = 2.0e3

        "--order"
            help = "Polynomial order for spectral elements (ignored if using config)"
            arg_type = Int
            default = 4

        "--smooth"
            help = "Number of mesh smoothing iterations (ignored if using config)"
            arg_type = Int
            default = 50
    end

    return parse_args(s)
end

function load_mesh_params_from_toml(config_file::String)
    """Load mesh parameters from TOML configuration file."""
    if !isfile(config_file)
        error("Config file not found: $config_file")
    end

    config = TOML.parsefile(config_file)

    # Extract parameters with error checking
    try
        params = MeshParameters(
            Lx = config["domain"]["Lx"],
            Ly = config["domain"]["Ly"],
            fault_start_y = config["fault"]["start_y"],
            background_grid_size = config["meshing"]["background_grid_size"],
            fault_refinement_h = config["meshing"]["fault_refinement_h"],
            fault_refinement_width = config["meshing"]["fault_refinement_width"],
            polynomial_order = config["meshing"]["polynomial_order"],
            smoothing_iterations = config["meshing"]["smoothing_iterations"]
        )
        return params, get(get(config, "output", Dict()), "mesh_name", "unstructured")
    catch e
        error("Error parsing config file $config_file: $e")
    end
end

function main()
    # Parse command line arguments
    args = parse_commandline()

    config_or_name = args["config_or_name"]

    # Determine if we're using a config file or direct parameters
    using_config = endswith(config_or_name, ".toml")

    if using_config
        # Config file mode
        if isnothing(args["simulation_name"])
            error("When using config file, you must provide simulation_name as second argument")
        end

        simulation_name = args["simulation_name"]
        config_file = config_or_name

        println("Loading mesh configuration from: $config_file")
        params, mesh_name = load_mesh_params_from_toml(config_file)

    else
        # Direct parameters mode
        simulation_name = config_or_name
        mesh_name = "unstructured"

        params = MeshParameters(
            Lx = args["Lx"],
            Ly = args["Ly"],
            fault_start_y = args["fault-y"],
            background_grid_size = args["background"],
            fault_refinement_h = args["refine"],
            fault_refinement_width = args["refine-width"],
            polynomial_order = args["order"],
            smoothing_iterations = args["smooth"]
        )
    end

    output_dir = joinpath("data", "mesh", simulation_name)

    # Print configuration
    println("━"^70)
    println("SEAS-SEME Mesh Generation")
    println("━"^70)
    println("Simulation name:  $simulation_name")
    println("Output directory: $output_dir")
    println()
    println("Mesh Configuration:")
    println("  Domain size:           $(params.Lx/1e3) × $(params.Ly/1e3) km")
    println("  Fault start:           y = $(params.fault_start_y/1e3) km")
    println("  Background elem size:  $(params.background_grid_size/1e3) km")
    println("  Fault refinement:      $(params.fault_refinement_h/1e3) km")
    println("  Refinement width:      $(params.fault_refinement_width/1e3) km")
    println("  Polynomial order:      $(params.polynomial_order)")
    println("  Smoothing iterations:  $(params.smoothing_iterations)")
    println("━"^70)
    println()

    # Generate mesh
    try
        mesh_file = generate_unstructured_mesh(params, output_dir, mesh_name=mesh_name)

        println()
        println("━"^70)
        println("✓ Mesh generation successful!")
        println("━"^70)
        println()
        println("Generated files:")
        println("  Mesh file:    $mesh_file")
        println("  Control file: $(joinpath(output_dir, "$(mesh_name).control"))")
        println("  Viz file:     $(joinpath(output_dir, "$(mesh_name).tec"))")
        println()
        println("Next steps:")
        println("  1. Visualize mesh in ParaView: $(joinpath(output_dir, "$(mesh_name).tec"))")
        println("  2. Reference in config.toml:")
        println("     [mesh]")
        println("     file = \"$mesh_file\"")
        println("     polynomial_degree = $(params.polynomial_order)")
        println("━"^70)

        return 0
    catch e
        println()
        println("━"^70)
        println("✗ Mesh generation failed!")
        println("━"^70)
        println()
        println("Error: $e")
        if isa(e, InterruptException)
            println("Interrupted by user")
        else
            println()
            println("Stack trace:")
            showerror(stdout, e, catch_backtrace())
        end
        println()
        return 1
    end
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
