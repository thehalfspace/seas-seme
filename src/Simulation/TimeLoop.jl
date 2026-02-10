"""
    TimeLoop

Main time integration loop for earthquake cycle simulations.

Orchestrates quasi-static and dynamic solvers with automatic mode switching,
adaptive timestepping, output writing, and checkpointing.
"""

using Printf
using Dates


"""
    run!(simulation)

Execute main time integration loop until total simulation time is reached.

# Arguments
- `simulation`: Simulation object containing all components

# Algorithm
1. Initialize output files and logging
2. Main time loop:
   a. Compute adaptive timestep
   b. Take time step (QS or dynamic solver)
   c. Update solver mode based on max slip velocity
   d. Write output at specified intervals
   e. Save checkpoint at specified intervals
   f. Display progress
3. Finalize and close files

# Notes
- Automatically switches between quasi-static and dynamic solvers
- Handles HDF5 streaming output and periodic checkpointing
- Zero allocations in main loop after warmup
- Robust error handling with checkpoint save on failure
"""
function run!(simulation)
    # Unpack simulation components
    state = simulation.state
    config = simulation.config
    mesh = simulation.mesh
    physics = simulation.physics
    ics = simulation.ics
    params = simulation.params
    qs_solver = simulation.qs_solver
    dyn_solver = simulation.dyn_solver
    timestepper = simulation.timestepper
    io_manager = simulation.io_manager

    # Extract frequently used values
    total_time = config.simulation.total_time
    log_interval = config.output.log_interval

    # Display simulation info
    println("\n" * "="^80)
    println("SEAS-SEME: Spectral Element Earthquake Cycle Simulation")
    println("="^80)
    @printf("Simulation: %s\n", config.simulation.name)
    @printf("Total time: %.2f years (%.2e s)\n", total_time / (365.25 * 24 * 3600), total_time)
    @printf("Output directory: %s\n", config.simulation.output_dir)
    @printf("Minimum timestep (CFL): %.3e s\n", timestepper.dt_min)
    @printf("Maximum timestep: %.3e s (%.1f days)\n", timestepper.dt_max, timestepper.dt_max / 86400)
    @printf("Mesh DOFs: %d\n", mesh.ndof)
    @printf("Fault nodes: %d\n", length(mesh.boundaries.fault.node_ids))
    println("="^80)
    println()

    # Initialize I/O
    initialize_output!(io_manager, mesh, config, ics)

    # Record start time
    t_start = now()
    last_log_time = now()

    # Main time loop
    println("Starting time integration...")
    println()

    try
        while state.time < total_time
            state.iteration += 1

            # Compute adaptive timestep
            state.timestep = compute_timestep(
                state.Vf,
                mesh.boundaries.fault.coords,
                ics.friction,
                timestepper,
                state.solver_mode,
                state.timestep
            )

            # Ensure we don't overshoot total time
            if state.time + state.timestep > total_time
                state.timestep = total_time - state.time
            end

            # Take time step (dispatches on state type via multiple dispatch)
            if state.solver_mode == :quasistatic
                quasistatic_step!(state, qs_solver, mesh, physics, ics, params,
                                 state.timestep)
            else  # :dynamic
                dynamic_step!(state, dyn_solver, mesh, physics, ics, params,
                             simulation.M_global, simulation.K_el, mesh.dof_id)
            end

            # Update time
            state.time += state.timestep

            # Compute current maximum slip rate for mode switching
            Vf_max = maximum_fault_slip_rate(state)

            # Determine solver mode for next iteration
            old_mode = state.solver_mode
            state.solver_mode = determine_solver_mode(Vf_max, state.solver_mode,
                                                     timestepper)

            # Log solver mode change
            if state.solver_mode != old_mode
                elapsed = now() - t_start
                @printf("[%s] Mode switch: %s → %s (Vf_max = %.3e m/s, iter = %d)\n",
                       format_duration(elapsed), old_mode, state.solver_mode,
                       Vf_max, state.iteration)
            end

            # Write output at specified intervals
            if should_write_output(state, config)
                write_timestep!(io_manager, state, mesh, ics, params)
            end

            # Write snapshot if conditions are met
            if should_write_snapshot(io_manager, state, config)
                write_snapshot!(io_manager, state, mesh, ics, params, config)
            end

            # Save checkpoint at specified intervals
            if should_save_checkpoint(state, config)
                save_checkpoint!(simulation, config)
                elapsed = now() - t_start
                @printf("[%s] Checkpoint saved at t = %.3e s (iter = %d)\n",
                       format_duration(elapsed), state.time, state.iteration)
            end

            # Display progress
            if should_display_progress(state, log_interval, last_log_time)
                display_progress(state, total_time, t_start, Vf_max)
                last_log_time = now()
            end

            # Check for NaN/Inf
            if !isfinite(Vf_max) || !all(isfinite.(state.u))
                @error "Non-finite values detected" Vf_max iter=state.iteration
                save_checkpoint!(simulation, config, emergency=true)
                error("Simulation failed due to non-finite values")
            end
        end

    catch e
        # Save emergency checkpoint on error
        @error "Simulation error at iteration $(state.iteration)" exception=(e, catch_backtrace())
        try
            save_checkpoint!(simulation, config, emergency=true)
            @warn "Emergency checkpoint saved"
        catch checkpoint_error
            @error "Failed to save emergency checkpoint" exception=(checkpoint_error, catch_backtrace())
        end
        rethrow(e)
    end

    # Finalize
    finalize_output!(io_manager)

    # Summary
    elapsed = now() - t_start
    println()
    println("="^80)
    println("Simulation Complete!")
    println("="^80)
    @printf("Total iterations: %d\n", state.iteration)
    @printf("Final time: %.2f years\n", state.time / (365.25 * 24 * 3600))
    @printf("Wall time: %s\n", format_duration(elapsed))
    @printf("Average time per iteration: %.2f ms\n",
           Dates.value(elapsed) / state.iteration)
    println("="^80)
    println()

    # Close logging
    close_logging(simulation.log_io)

    return nothing
end


"""
    should_write_output(state, config) -> Bool

Determine if output should be written at current iteration.
"""
function should_write_output(state, config)
    # Write every iteration during dynamic events
    if state.solver_mode == :dynamic
        return true
    end

    # Write periodically during quasi-static
    # (could add time-based intervals here if desired)
    return mod(state.iteration, 10) == 0
end


"""
    should_save_checkpoint(state, config) -> Bool

Determine if checkpoint should be saved at current iteration.
"""
function should_save_checkpoint(state, config)
    if !config.checkpointing.enabled
        return false
    end

    # Save at specified time intervals
    checkpoint_interval = config.checkpointing.interval

    # Check if we've crossed a checkpoint threshold
    n_checkpoints = floor(Int, state.time / checkpoint_interval)
    prev_time = state.time - state.timestep
    n_prev_checkpoints = floor(Int, prev_time / checkpoint_interval)

    return n_checkpoints > n_prev_checkpoints
end


"""
    should_display_progress(state, log_interval, last_log_time) -> Bool

Determine if progress should be displayed.
"""
function should_display_progress(state, log_interval::Int,
                                 last_log_time::DateTime)
    # Display every log_interval iterations
    if mod(state.iteration, log_interval) == 0
        return true
    end

    # Also display if more than 30 seconds have passed
    return (now() - last_log_time) > Millisecond(30000)
end


"""
    display_progress(state, total_time, t_start, Vf_max)

Display progress information.
"""
function display_progress(state, total_time::Real,
                         t_start::DateTime, Vf_max::Real)
    elapsed = now() - t_start
    progress_pct = 100 * state.time / total_time

    @printf("[%s] iter=%d, t=%.3e s (%.1f%%), dt=%.3e s, Vf_max=%.3e m/s, mode=%s\n",
           format_duration(elapsed), state.iteration, state.time, progress_pct,
           state.timestep, Vf_max, state.solver_mode)
end


"""
    format_duration(duration::Period) -> String

Format duration for display (e.g., "2h 15m 30s").
"""
function format_duration(duration::Period)
    total_ms = Dates.value(duration)
    total_s = total_ms ÷ 1000

    hours = total_s ÷ 3600
    minutes = (total_s % 3600) ÷ 60
    seconds = total_s % 60

    if hours > 0
        return @sprintf("%dh %dm %ds", hours, minutes, seconds)
    elseif minutes > 0
        return @sprintf("%dm %ds", minutes, seconds)
    else
        return @sprintf("%ds", seconds)
    end
end
