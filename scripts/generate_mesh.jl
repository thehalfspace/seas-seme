"""
Mesh Generation Script for SEAS-SEME

This script generates unstructured quadrilateral meshes for earthquake cycle
simulations using HOHQMesh.jl.

# Usage

## From configuration file (recommended):
```bash
julia scripts/generate_mesh.jl config/mesh_default.toml
```

## From command line with custom parameters:
```bash
julia scripts/generate_mesh.jl simulation_name --Lx 60e3 --fault-y 30e3 --refine 0.1e3
```

# Command Line Arguments

Positional:
- `config_or_name`: Either a .toml config file path OR simulation name

When using config file:
- `simulation_name`: Optional display name (if not provided, uses directory name from config)

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
julia scripts/generate_mesh.jl config/mesh_default.toml
```

## Using command line (for quick tests):
```bash
julia scripts/generate_mesh.jl test_mesh --refine 0.1e3
```

# Output

Creates mesh files in the directory specified by `output_dir` in the config:
- `{mesh_name}.control`: HOHQMesh control file
- `{mesh_name}.mesh`: Binary mesh file (use in config.toml)
- `{mesh_name}.tec`: Tecplot file (visualize in ParaView)

The mesh file path should be referenced in your simulation config.toml:
```toml
[mesh]
file = "data/mesh/homogeneous/unstructured_02.mesh"
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
            julia scripts/generate_mesh.jl config/mesh_default.toml
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

        output_config = get(config, "output", Dict())
        mesh_name = get(output_config, "mesh_name", "unstructured")
        output_dir = get(output_config, "output_dir", "data/mesh")

        return params, mesh_name, output_dir
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
        config_file = config_or_name

        println("Loading mesh configuration from: $config_file")
        params, mesh_name, output_dir = load_mesh_params_from_toml(config_file)

        # Simulation name is optional now (used only for display)
        simulation_name = if !isnothing(args["simulation_name"])
            args["simulation_name"]
        else
            basename(output_dir)  # Use directory name as simulation name
        end

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

        # Use default output directory for direct mode
        output_dir = joinpath("data", "mesh", simulation_name)
    end

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
