"""
    Checkpointing

Save and load simulation checkpoints using JLD2.

Features:
- Complete simulation state serialization
- Automatic checkpoint rotation (keep last N)
- Emergency checkpoints on errors
- Fast restart from checkpoint
"""

using JLD2
using Dates
using Printf


"""
    save_checkpoint!(simulation, config; emergency=false)

Save simulation checkpoint to disk.

# Arguments
- `simulation`: Simulation object to save
- `config::SimulationConfig`: Configuration with checkpoint settings
- `emergency::Bool`: Mark as emergency checkpoint (default: false)

# Saves
- Complete simulation state
- Solver objects (including AMG preconditioner)
- Timestamp and metadata

# Checkpoint Rotation
Automatically removes old checkpoints, keeping only the last N
(specified by `config.checkpointing.keep_last`).

# Filename Format
- Normal: `checkpoint_iter_<iteration>.jld2`
- Emergency: `checkpoint_emergency_iter_<iteration>.jld2`
"""
function save_checkpoint!(simulation, config; emergency::Bool=false)
    if !config.checkpointing.enabled && !emergency
        return
    end

    # Create checkpoint directory
    ckpt_dir = config.checkpointing.directory
    mkpath(ckpt_dir)

    # Generate filename
    iter = simulation.state.iteration
    if emergency
        filename = "checkpoint_emergency_iter_$(iter).jld2"
    else
        filename = "checkpoint_iter_$(iter).jld2"
    end
    filepath = joinpath(ckpt_dir, filename)

    # Save checkpoint
    jldsave(filepath;
        # Configuration
        config = simulation.config,

        # Simulation state
        state = simulation.state,

        # Mesh and geometry (lightweight - just indices and coordinates)
        mesh_ndof = simulation.mesh.ndof,
        mesh_n_elements = simulation.mesh.n_elements,
        mesh_polynomial_degree = simulation.mesh.polynomial_degree,

        # Initial conditions and parameters
        ics = simulation.ics,
        params = simulation.params,

        # Solvers (including AMG preconditioner - expensive to rebuild!)
        qs_solver = simulation.qs_solver,
        dyn_solver = simulation.dyn_solver,
        timestepper = simulation.timestepper,

        # Matrices
        M_global = simulation.M_global,
        K_el = simulation.K_el,

        # Metadata
        save_time = now(),
        julia_version = string(VERSION),
        emergency = emergency
    )

    # Cleanup old checkpoints (only for normal checkpoints)
    if !emergency
        cleanup_old_checkpoints!(ckpt_dir, config.checkpointing.keep_last)
    end

    return filepath
end


"""
    load_checkpoint(filepath::String, mesh) -> Simulation

Load simulation from checkpoint file.

# Arguments
- `filepath::String`: Path to checkpoint file
- `mesh`: Reconstructed mesh object

# Returns
- `Simulation`: Loaded simulation ready to resume

# Notes
The mesh must be reconstructed separately (from config) because it
contains references to external mesh files. All other components
are fully serialized in the checkpoint.

# Example
```julia
config = load_config("config.toml")
mesh = build_mesh(config.mesh)
sim = load_checkpoint("checkpoint.jld2", mesh)
run!(sim)  # Resume from checkpoint
```
"""
function load_checkpoint(filepath::String, mesh)
    @info "Loading checkpoint" filepath=filepath

    # Load checkpoint data
    ckpt = load(filepath)

    # Extract components
    config = ckpt["config"]
    state = ckpt["state"]
    ics = ckpt["ics"]
    params = ckpt["params"]
    qs_solver = ckpt["qs_solver"]
    dyn_solver = ckpt["dyn_solver"]
    timestepper = ckpt["timestepper"]
    M_global = ckpt["M_global"]
    K_el = ckpt["K_el"]

    # Verify mesh compatibility
    if mesh.ndof != ckpt["mesh_ndof"] ||
       mesh.n_elements != ckpt["mesh_n_elements"] ||
       mesh.polynomial_degree != ckpt["mesh_polynomial_degree"]
        error("Mesh mismatch: checkpoint mesh incompatible with current mesh")
    end

    # Reconstruct physics
    physics = MaterialProperties(
        config.physics.density,
        config.physics.shear_velocity
    )

    # Create I/O manager (will be reinitialized)
    io_manager = create_io_manager(config, mesh)

    # Reconstruct simulation
    T = eltype(state.u)
    simulation = Simulation{T}(
        config, mesh, physics, ics, params,
        state, qs_solver, dyn_solver, timestepper,
        io_manager, M_global, K_el
    )

    # Display checkpoint info
    save_time = ckpt["save_time"]
    emergency = get(ckpt, "emergency", false)

    @info "Checkpoint loaded successfully" emergency=emergency save_time=save_time
    @printf("  Resuming from: t = %.3e s, iteration = %d\n",
           state.time, state.iteration)
    @printf("  Solver mode: %s\n", state.solver_mode)
    @printf("  Timestep: %.3e s\n", state.timestep)

    return simulation
end


"""
    cleanup_old_checkpoints!(directory::String, keep_last::Int)

Remove old checkpoints, keeping only the most recent N.

# Arguments
- `directory::String`: Checkpoint directory
- `keep_last::Int`: Number of checkpoints to retain

# Notes
- Only removes normal checkpoints (not emergency checkpoints)
- Sorts by modification time (most recent kept)
"""
function cleanup_old_checkpoints!(directory::String, keep_last::Int)
    # Find all normal checkpoint files
    all_files = readdir(directory, join=true)
    checkpoint_files = filter(f -> startswith(basename(f), "checkpoint_iter_") &&
                                  endswith(f, ".jld2"), all_files)

    # Sort by modification time (newest first)
    sort!(checkpoint_files, by=mtime, rev=true)

    # Remove old checkpoints
    if length(checkpoint_files) > keep_last
        for f in checkpoint_files[keep_last+1:end]
            rm(f)
            @debug "Removed old checkpoint" file=basename(f)
        end
    end

    return nothing
end


"""
    find_latest_checkpoint(directory::String) -> Union{String, Nothing}

Find the most recent checkpoint file in directory.

# Arguments
- `directory::String`: Checkpoint directory

# Returns
- `String`: Path to latest checkpoint, or `nothing` if no checkpoints found

# Notes
Searches for normal checkpoints only (not emergency checkpoints).
"""
function find_latest_checkpoint(directory::String)
    if !isdir(directory)
        return nothing
    end

    # Find all normal checkpoint files
    all_files = readdir(directory, join=true)
    checkpoint_files = filter(f -> startswith(basename(f), "checkpoint_iter_") &&
                                  endswith(f, ".jld2"), all_files)

    if isempty(checkpoint_files)
        return nothing
    end

    # Return most recent
    return checkpoint_files[argmax(mtime.(checkpoint_files))]
end
