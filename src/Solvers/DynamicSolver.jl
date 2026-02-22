"""
    DynamicSolver

Explicit dynamic earthquake solver using leap-frog time integration.

Implements the algorithm from Kaneko et al. (2008) Appendix B for wave propagation
during rapid earthquake slip. Uses Newton-Raphson search to enforce rate-state
friction boundary conditions on the fault.
"""

# Note: All dependencies are imported at the main module level (SEAS_SEME.jl).
# Functions from other submodules are available in the same namespace.


"""
    DynamicSolver{T}

Explicit dynamic solver configuration.

# Fields
- `dt_min::T`: Minimum timestep from CFL condition [s]
- `verbose::Bool`: Print diagnostic information
"""
struct DynamicSolver{T<:AbstractFloat}
    dt_min::T
    verbose::Bool
end


"""
    DynamicSolver(dt_min; verbose=false)

Construct dynamic solver with specified minimum timestep.

# Arguments
- `dt_min`: CFL timestep [s] (typically ≈ dx / (√2 * vs) for stability)
- `verbose`: Print diagnostics (default: false)
"""
function DynamicSolver(dt_min::T; verbose::Bool=false) where T<:AbstractFloat
    return DynamicSolver{T}(dt_min, verbose)
end


"""
    dynamic_step!(state, solver, mesh, physics, ics, params, M_global, K_el, dof_id)

Perform one dynamic time step using leap-frog integration.

# Algorithm (Kaneko et al. 2008, Appendix B)
1. Update displacement: u[n+1] = u[n] + dt*v[n] + 0.5*dt²*a[n]
2. Partial velocity update: v[n+1/2] = v[n] + 0.5*dt*a[n]
3. Compute internal forces: f_int = -K*u[n+1]
4. Apply absorbing boundary conditions
5. Compute stick traction (free velocity) on fault
6. Update fault state variable
7. Newton-Raphson search to satisfy:
   - Stick traction equation
   - Rate-state friction law
   (Two NR iterations for accuracy)
8. Apply fault traction to force vector
9. Compute acceleration: a[n+1] = f / M
10. Complete velocity update: v[n+1] = v[n+1/2] + 0.5*dt*a[n+1]

# Arguments
- `state`: SimulationState (modified in place)
- `solver`: DynamicSolver
- `mesh`: UnstructuredSEMesh
- `physics`: Material properties
- `ics`: Initial conditions
- `params`: Simulation parameters
- `M_global`: Global mass matrix (diagonal)
- `K_el`: Elemental stiffness matrices
- `dof_id`: DOF connectivity

# Modifies
- `state.u`, `state.v`, `state.a`: Displacement, velocity, acceleration
- `state.τf`, `state.ψ`, `state.Vf`: Fault stress, state, slip rate
- `state.fault_vfree`: Stick traction workspace
"""
function dynamic_step!(state, solver::DynamicSolver{T}, mesh, physics, ics,
                      params, M_global::Vector{T}, K_el::Array{T,3},
                      dof_id::Array{Int,3}) where T<:AbstractFloat

    dt = solver.dt_min

    # Extract boundary info
    fault_id = mesh.boundaries.fault.node_ids
    creep_id = mesh.boundaries.creep.node_ids
    absorbing_id = mesh.boundaries.absorbing.node_ids
    fault_matrix = mesh.boundaries.fault.matrix

    # Compute fault impedance: Z = M / (fault_matrix * dt)
    # See Kaneko et al. (2008) for derivation
    fault_z = M_global[fault_id] ./ (fault_matrix .* dt)

    absorb_matrix = mesh.boundaries.absorbing.matrix

    # Store previous state
    copyto!(state.u_prev, state.u)
    copyto!(state.v_prev, state.v)

    # Step 1: Update displacement using leap-frog
    # u[n+1] = u[n] + dt*v[n] + 0.5*dt²*a[n]
    state.u .= state.u .+ dt .* state.v .+ (0.5 * dt^2) .* state.a

    # Step 2: Partial velocity update
    # v[n+1/2] = v[n] + 0.5*dt*a[n]
    state.v .= state.v .+ (0.5 * dt) .* state.a
    fill!(state.a, zero(T))

    # Step 3: Internal forces -K*u[n+1]
    state.a .= apply_stiffness!(state.a, state.u, K_el, dof_id, mesh.n_elements)
    state.a .= -state.a

    # Enforce zero force on creep boundary
    state.a[creep_id] .= 0

    # Apply absorbing boundary conditions (Lysmer dampers)
    state.a[absorbing_id] .= state.a[absorbing_id] .-
                             (absorb_matrix .* state.v[absorbing_id])

    #--------------------------------------------------
    # Rate-state fault boundary
    #--------------------------------------------------

    # Step 4: Compute stick traction (free velocity if no friction)
    # FaultVFree = 2*v[n+1/2] + dt*a_internal/M
    state.fault_vfree .= 2 .* state.v[fault_id] .+
                        dt .* state.a[fault_id] ./ M_global[fault_id]

    # Fault slip rate from previous timestep (for state update)
    Vf_prev = 2 .* state.v_prev[fault_id] .+ params.Vpl

    # Workspace for two-iteration NR search
    Vf_first_iter = similar(state.Vf)

    # Step 5-7: Solve nonlinear fault boundary equations
    for i in eachindex(fault_id)
        # Save ψ before update so we can restore on NR failure
        ψ_saved = state.ψ[i]

        # Update state variable using slip rate from previous step
        state.ψ[i] = state_time_evolution(state.ψ[i], Vf_prev[i], dt,
                                         ics.friction.Lc[i], params.Vo)

        # Store slip rate for second iteration
        Vf_first_iter[i] = Vf_prev[i]

        # CRITICAL FIX: Convert perturbation stress → total stress for nr_search
        # state.τf[i] stores τ_perturbation = τ_total - τ_init from previous timestep
        # But nr_search expects τ_total as initial guess
        τ_total_guess = state.τf[i] + ics.τo[i]

        # Newton-Raphson search (1st iteration)
        state.Vf[i], state.τf[i] = nr_search(
            τ_total_guess,        # FIXED: Use total stress as initial guess
            params.fo,
            params.Vo,
            ics.friction.a[i],
            ics.friction.b[i],
            ics.σo[i],
            ics.τo[i],
            state.ψ[i],
            fault_z[i],
            state.fault_vfree[i]
        )

        # Check for NaN (indicates NR failure) — restore ψ and continue
        if isnan(state.Vf[i]) || isnan(state.τf[i])
            state.ψ[i] = ψ_saved
            state.Vf[i] = Vf_prev[i]
            state.τf[i] = τ_total_guess - ics.τo[i]
            @warn "NR search produced NaN, restoring previous state" location=i iter=state.iteration
            continue
        end

        # 2nd iteration: Update state with average slip rate for accuracy
        # (See Kaneko et al. 2008 - single iteration is less accurate)
        ψ_after_first = state.ψ[i]
        avg_slip_rate = 0.5 * (state.Vf[i] + Vf_first_iter[i])
        state.ψ[i] = state_time_evolution(state.ψ[i], avg_slip_rate, dt,
                                         ics.friction.Lc[i], params.Vo)

        # Newton-Raphson search (2nd iteration) — capture locally to guard NaN
        # state.τf[i] already contains total stress from 1st iteration, no conversion needed
        Vf_2nd, τf_2nd = nr_search(
            state.τf[i],          # Already total stress from 1st iteration
            params.fo,
            params.Vo,
            ics.friction.a[i],
            ics.friction.b[i],
            ics.σo[i],
            ics.τo[i],
            state.ψ[i],
            fault_z[i],
            state.fault_vfree[i]
        )

        if isnan(Vf_2nd) || isnan(τf_2nd)
            # 2nd NR failed — keep 1st iteration results and restore ψ
            state.ψ[i] = ψ_after_first
            @warn "NR 2nd iteration produced NaN, keeping 1st iteration result" location=i
        else
            state.Vf[i] = Vf_2nd
            state.τf[i] = τf_2nd
        end
    end

    # Convert fault stress to perturbation from initial stress
    # (NR search works with total stress = τ_init + τ_perturbation)
    state.τf .= state.τf .- ics.τo

    # Step 8: Apply fault traction to force vector
    state.a[fault_id] .= state.a[fault_id] .- fault_matrix .* state.τf

    #--------------------------------------------------
    # End of fault boundary
    #--------------------------------------------------

    # Step 9: Compute acceleration from forces
    # a[n+1] = f / M
    state.a .= state.a ./ M_global

    # Step 10: Complete velocity update
    # v[n+1] = v[n+1/2] + 0.5*dt*a[n+1]
    state.v .= state.v .+ (0.5 * dt) .* state.a

    # Enforce zero velocity on creep boundary
    state.v[creep_id] .= 0
    state.a[creep_id] .= 0

    return nothing
end
