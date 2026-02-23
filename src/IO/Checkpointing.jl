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
    _state_to_cpu(state)

Download GPU state arrays to CPU for serialization. Returns a CPU-only copy of state.
For CPU states, returns the state unchanged.
"""
function _state_to_cpu(state::SimulationState{T, V}) where {T, V}
    if V <: Vector
        return state   # already CPU
    end
    # GPU state: download all array fields
    return SimulationState{T, Vector{T}}(
        Array(state.u), Array(state.v), Array(state.a),
        Array(state.τf), Array(state.ψ), Array(state.Vf), Array(state.cum_slip),
        Array(state.u_prev), Array(state.v_prev), Array(state.f), Array(state.fault_vfree),
        state.time, state.timestep, state.iteration, state.solver_mode
    )
end

function _state_to_cpu(state::SimulationStatePlaneStrain{T, V}) where {T, V}
    if V <: Vector
        return state
    end
    return SimulationStatePlaneStrain{T, Vector{T}}(
        Array(state.u), Array(state.v), Array(state.a),
        Array(state.τf), Array(state.ψ), Array(state.Vf), Array(state.cum_slip),
        Array(state.σn_perturbation),
        Array(state.u_prev), Array(state.v_prev), Array(state.f), Array(state.fault_vfree),
        state.fault_tangent, state.fault_normal,
        state.time, state.timestep, state.iteration, state.solver_mode,
        state.ndof, state.nfault
    )
end


"""
    save_checkpoint!(simulation, config; emergency=false)

Save simulation checkpoint to disk.

# Arguments
- `simulation`: Simulation object to save
- `config::SimulationConfig`: Configuration with checkpoint settings
- `emergency::Bool`: Mark as emergency checkpoint (default: false)

# Saves
- Complete simulation state (downloaded to CPU if on GPU)
- Solver objects (CPU-side AMG or CPU copy of GPU solver config)
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

    # Download GPU arrays to CPU before serialization
    state_cpu    = _state_to_cpu(simulation.state)
    M_global_cpu = Array(simulation.M_global)
    weights_cpu  = cpu(simulation.weights)
    H_cpu        = Matrix(simulation.H)
    Ht_cpu       = Matrix(simulation.Ht)

    # Note: GPU AMGX solvers are not JLD2-serializable.
    # Save a marker so load_checkpoint knows to rebuild on GPU.
    qs_solver_save = simulation.use_gpu ? nothing : simulation.qs_solver

    # Save checkpoint
    jldsave(filepath;
        # Configuration
        config = simulation.config,

        # Simulation state (CPU)
        state = state_cpu,

        # Mesh and geometry (lightweight - just indices and coordinates)
        mesh_ndof = simulation.mesh.ndof,
        mesh_n_elements = simulation.mesh.n_elements,
        mesh_polynomial_degree = simulation.mesh.polynomial_degree,

        # Initial conditions and parameters
        ics = simulation.ics,
        params = simulation.params,

        # Solvers (AMGX GPU solvers are rebuilt on load; CPU AMG is serialized)
        qs_solver = qs_solver_save,
        dyn_solver = simulation.dyn_solver,
        timestepper = simulation.timestepper,

        # Matrices (CPU)
        M_global = M_global_cpu,
        weights  = weights_cpu,
        H        = H_cpu,
        Ht       = Ht_cpu,

        # GPU flag (needed to know how to restore)
        use_gpu = simulation.use_gpu,

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
    config     = ckpt["config"]
    state      = ckpt["state"]
    ics        = ckpt["ics"]
    params     = ckpt["params"]
    qs_solver  = ckpt["qs_solver"]
    dyn_solver = ckpt["dyn_solver"]
    timestepper = ckpt["timestepper"]
    M_global   = ckpt["M_global"]
    weights    = ckpt["weights"]
    H          = ckpt["H"]
    Ht         = ckpt["Ht"]
    use_gpu    = get(ckpt, "use_gpu", false)

    # Verify mesh compatibility
    if mesh.ndof != ckpt["mesh_ndof"] ||
       mesh.n_elements != ckpt["mesh_n_elements"] ||
       mesh.polynomial_degree != ckpt["mesh_polynomial_degree"]
        error("Mesh mismatch: checkpoint mesh incompatible with current mesh")
    end

    # Reconstruct physics (formulation-aware)
    if config.physics.formulation == :plane_strain
        physics = MaterialProperties(
            config.physics.density,
            config.physics.shear_velocity,
            config.physics.poisson_ratio
        )
    else
        physics = MaterialProperties(
            config.physics.density,
            config.physics.shear_velocity
        )
    end

    # GPU: re-upload arrays and rebuild AMGX solver
    dof_id = mesh.dof_id
    if use_gpu
        @info "Restoring GPU arrays from checkpoint..."
        M_global = CuArray(M_global)
        weights  = gpu(weights)
        H        = CuMatrix(H)
        Ht       = CuMatrix(Ht)
        dof_id   = CuArray(dof_id)

        # Rebuild AMGX solver (cannot be serialized)
        fltni = _reconstruct_fltni(mesh, config)
        if config.physics.formulation == :plane_strain
            w_cpu = cpu(weights)  # temporarily cpu for assembly
            qs_solver = build_quasistatic_solver_gpu_plane_strain(
                w_cpu, Matrix(H), Matrix(Ht), mesh.dof_id, mesh, fltni
            )
        else
            w_cpu = cpu(weights)
            qs_solver = build_quasistatic_solver_gpu(
                w_cpu, Matrix(H), Matrix(Ht), mesh.dof_id, mesh, fltni
            )
        end

        # Re-upload state arrays to GPU
        state = _state_to_gpu(state)
    end

    # Create I/O manager (will be reinitialized)
    io_manager = create_io_manager(config, mesh)

    # Setup logging for resumed simulation
    _, log_io = setup_logging(config)

    # Reconstruct simulation
    simulation = Simulation(
        config, mesh, physics, ics, params,
        state, qs_solver, dyn_solver, timestepper,
        io_manager, log_io, M_global, weights, H, Ht, dof_id, use_gpu
    )

    # Display checkpoint info
    save_time = ckpt["save_time"]
    emergency = get(ckpt, "emergency", false)

    @info "Checkpoint loaded successfully" emergency=emergency save_time=save_time use_gpu=use_gpu
    @printf("  Resuming from: t = %.3e s, iteration = %d\n",
           state.time, state.iteration)
    @printf("  Solver mode: %s\n", state.solver_mode)
    @printf("  Timestep: %.3e s\n", state.timestep)

    return simulation
end


"""
    _state_to_gpu(state)

Upload CPU state arrays to GPU after loading from checkpoint.
"""
function _state_to_gpu(state::SimulationState{T, Vector{T}}) where T
    return SimulationState{T, CuVector{T}}(
        CuArray(state.u), CuArray(state.v), CuArray(state.a),
        CuArray(state.τf), CuArray(state.ψ), CuArray(state.Vf), CuArray(state.cum_slip),
        CuArray(state.u_prev), CuArray(state.v_prev),
        CuArray(state.f), CuArray(state.fault_vfree),
        state.time, state.timestep, state.iteration, state.solver_mode
    )
end

function _state_to_gpu(state::SimulationStatePlaneStrain{T, Vector{T}}) where T
    return SimulationStatePlaneStrain{T, CuVector{T}}(
        CuArray(state.u), CuArray(state.v), CuArray(state.a),
        CuArray(state.τf), CuArray(state.ψ), CuArray(state.Vf), CuArray(state.cum_slip),
        CuArray(state.σn_perturbation),
        CuArray(state.u_prev), CuArray(state.v_prev),
        CuArray(state.f), CuArray(state.fault_vfree),
        state.fault_tangent, state.fault_normal,
        state.time, state.timestep, state.iteration, state.solver_mode,
        state.ndof, state.nfault
    )
end

# No-op: already GPU or unknown type
_state_to_gpu(state) = state


"""
    _reconstruct_fltni(mesh, config)

Reconstruct the free DOF index vector from mesh and config.
Mirrors the logic in build_simulation_antiplane/plane_strain.
"""
function _reconstruct_fltni(mesh, config)
    ndof = mesh.ndof
    interface_id_spatial = sort(unique(vcat(mesh.boundaries.creep.node_ids,
                                            mesh.boundaries.fault.node_ids)))
    if config.physics.formulation == :plane_strain
        interface_id_2n = vcat(interface_id_spatial, ndof .+ interface_id_spatial)
        fltni = collect(1:2*ndof)
        deleteat!(fltni, sort(interface_id_2n))
    else
        fltni = collect(1:ndof)
        deleteat!(fltni, sort(unique(interface_id_spatial)))
    end
    return fltni
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
