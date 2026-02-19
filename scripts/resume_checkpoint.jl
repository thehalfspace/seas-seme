#!/usr/bin/env julia
"""
resume_checkpoint.jl - Resume a simulation from a checkpoint file.

# Usage
    julia --project=. scripts/resume_checkpoint.jl <config_file> <checkpoint_file>

# Example
    julia --project=. scripts/resume_checkpoint.jl config/dip_slip_2d.toml \\
        data/dip_slip_2d_test/checkpoints/checkpoint_emergency_iter_187198.jld2
"""

using SEAS_SEME
using Printf

function main()
    if length(ARGS) < 2
        println("Usage: julia resume_checkpoint.jl <config_file> <checkpoint_file>")
        exit(1)
    end

    config_file = ARGS[1]
    checkpoint_file = ARGS[2]

    isfile(config_file)     || (@error "Config not found" config_file;         exit(1))
    isfile(checkpoint_file) || (@error "Checkpoint not found" checkpoint_file; exit(1))

    println()
    println("="^80)
    println("SEAS-SEME: Resume from Checkpoint")
    println("="^80)
    @printf("Configuration: %s\n", config_file)
    @printf("Checkpoint:    %s\n", checkpoint_file)
    @printf("Julia threads: %d\n", Threads.nthreads())
    println("="^80)
    println()

    println("Loading configuration...")
    config = load_config(config_file)

    # build_simulation constructs the mesh (needed by load_checkpoint for compatibility check)
    println("Building simulation (for mesh)...")
    simulation = build_simulation(config)
    println("✓ Simulation built")

    println("Loading checkpoint...")
    simulation = load_checkpoint(checkpoint_file, simulation.mesh)
    println("✓ Checkpoint loaded")

    println("\nResuming simulation...")
    run!(simulation)

    println("\n✓ Simulation complete!")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
