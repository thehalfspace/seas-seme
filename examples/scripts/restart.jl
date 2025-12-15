#!/usr/bin/env julia

"""
restart.jl - Restart SEAS-SEME simulation from checkpoint

Resume a simulation from a saved checkpoint file.

# Usage

```bash
# Restart from specific checkpoint
julia --project=. examples/scripts/restart.jl <checkpoint_file>

# Restart from latest checkpoint in directory
julia --project=. examples/scripts/restart.jl --latest <checkpoint_dir>
```

# Example

```bash
# Restart from specific checkpoint
julia --project=. examples/scripts/restart.jl checkpoints/checkpoint_iter_50000.jld2

# Restart from latest checkpoint
julia --project=. examples/scripts/restart.jl --latest checkpoints/
```

# Notes

- The original mesh file must be accessible (specified in config)
- Output will be appended to existing HDF5 file
- New checkpoints will be saved as simulation continues
"""

using SEAS_SEME
using Printf


function main()
    # Parse command line arguments
    if length(ARGS) < 1
        println("Usage:")
        println("  julia restart.jl <checkpoint_file>")
        println("  julia restart.jl --latest <checkpoint_dir>")
        println()
        println("Example:")
        println("  julia restart.jl checkpoints/checkpoint_iter_50000.jld2")
        println("  julia restart.jl --latest checkpoints/")
        exit(1)
    end

    # Determine checkpoint file
    if ARGS[1] == "--latest"
        if length(ARGS) < 2
            @error "Must specify checkpoint directory with --latest"
            exit(1)
        end

        checkpoint_dir = ARGS[2]
        checkpoint_file = find_latest_checkpoint(checkpoint_dir)

        if checkpoint_file === nothing
            @error "No checkpoints found in directory" directory=checkpoint_dir
            exit(1)
        end

        println("Found latest checkpoint: ", basename(checkpoint_file))
    else
        checkpoint_file = ARGS[1]
    end

    # Verify checkpoint file exists
    if !isfile(checkpoint_file)
        @error "Checkpoint file not found" checkpoint_file
        exit(1)
    end

    println()
    println("="^80)
    println("SEAS-SEME: Restart from Checkpoint")
    println("="^80)
    @printf("Checkpoint: %s\n", checkpoint_file)
    @printf("Julia threads: %d\n", Threads.nthreads())
    println("="^80)
    println()

    # Load checkpoint (partial load to get config)
    println("Loading checkpoint configuration...")
    ckpt = JLD2.load(checkpoint_file)
    config = ckpt["config"]

    # Rebuild mesh (required for restart)
    println("\nRebuilding mesh...")
    mesh = build_mesh(config.mesh)
    println("✓ Mesh rebuilt")

    # Load full simulation from checkpoint
    println("\nLoading simulation state...")
    simulation = load_checkpoint(checkpoint_file, mesh)
    println("✓ Simulation loaded")

    # Resume simulation
    println("\nResuming simulation...")
    run!(simulation)

    println("\n✓ Simulation complete!")
    println()

    return 0
end


# Run main function
if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
