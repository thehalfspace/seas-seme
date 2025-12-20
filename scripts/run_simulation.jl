#!/usr/bin/env julia

"""
run_simulation.jl - Main entry point for SEAS-SEME simulations

Run a complete earthquake cycle simulation from a configuration file.

# Usage

```bash
julia --project=. examples/scripts/run_simulation.jl <config_file>
```

# Example

```bash
julia --project=. examples/scripts/run_simulation.jl examples/config/strike_slip_2d.toml
```

# Options

Set number of threads before running:
```bash
export JULIA_NUM_THREADS=16
julia --project=. examples/scripts/run_simulation.jl config.toml
```

# Output

- HDF5 file: `<output_dir>/<simulation_name>.h5`
- Checkpoints: `<checkpoint_dir>/checkpoint_iter_*.jld2`
- Log file: `<output_dir>/<simulation_name>.log`
"""

using SEAS_SEME
using Printf


function main()
    # Parse command line arguments
    if length(ARGS) < 1
        println("Usage: julia run_simulation.jl <config_file>")
        println()
        println("Example:")
        println("  julia run_simulation.jl examples/config/strike_slip_2d.toml")
        exit(1)
    end

    config_file = ARGS[1]

    # Verify config file exists
    if !isfile(config_file)
        @error "Configuration file not found" config_file
        exit(1)
    end

    println()
    println("="^80)
    println("SEAS-SEME: Spectral Element Method for Earthquake Cycle Simulations")
    println("="^80)
    @printf("Configuration: %s\n", config_file)
    @printf("Julia threads: %d\n", Threads.nthreads())
    println("="^80)
    println()

    # Load configuration
    println("Loading configuration...")
    config = load_config(config_file)
    println("✓ Configuration loaded")

    # Setup simulation directories and save config
    setup_simulation_directories(config, config_file)

    # Build simulation
    println("\nBuilding simulation...")
    simulation = build_simulation(config)
    println("✓ Simulation built")

    # Run simulation
    println("\nStarting simulation...")
    run!(simulation)

    println("\n✓ Simulation complete!")
    println()

    return 0
end


# Run main function
if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
