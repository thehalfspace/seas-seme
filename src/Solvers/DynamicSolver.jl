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
    dynamic_step!(state, solver, mesh, physics, ics, params, M_global, weights, H, Ht, dof_id)

Perform one dynamic time step using leap-frog integration.

# Algorithm (Kaneko et al. 2008, Appendix B)
1. Update displacement: u[n+1] = u[n] + dt*v[n] + 0.5*dt²*a[n]
2. Partial velocity update: v[n+1/2] = v[n] + 0.5*dt*a[n]
3. Compute internal forces: f_int = -K*u[n+1] (via tensor-product matvec)
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
- `weights::MetricWeightsAntiplane{T}`: Compact metric weights
- `H::Matrix{T}`: Derivative matrix
- `Ht::Matrix{T}`: H' (pre-transposed)
- `dof_id`: DOF connectivity

# Modifies
- `state.u`, `state.v`, `state.a`: Displacement, velocity, acceleration
- `state.τf`, `state.ψ`, `state.Vf`: Fault stress, state, slip rate
- `state.fault_vfree`: Stick traction workspace
"""
function dynamic_step!(state, solver::DynamicSolver{T}, mesh, physics, ics,
                      params, M_global::AbstractVector{T},
                      weights::MetricWeightsAntiplane{T},
                      H::AbstractMatrix{T}, Ht::AbstractMatrix{T},
                      dof_id::AbstractArray{<:Integer,3}) where T<:AbstractFloat

    dt = solver.dt_min

    # Extract boundary info
    fault_id = mesh.boundaries.fault.node_ids
    creep_id = mesh.boundaries.creep.node_ids
    absorbing_id = mesh.boundaries.absorbing.node_ids
    fault_matrix = mesh.boundaries.fault.matrix

    # Compute fault impedance: Z = M / (fault_matrix * dt)
    # See Kaneko et al. (2008) for derivation
    # Always keep fault_z on CPU (small array, needed for per-node NR loop)
    fault_z = Array(M_global[fault_id]) ./ (fault_matrix .* dt)

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

    # Step 3: Internal forces -K*u[n+1] via tensor-product matvec
    apply_stiffness!(state.a, state.u, weights, H, Ht, dof_id, mesh.n_elements)
    state.a .= -state.a

    # Enforce zero force on creep boundary
    state.a[creep_id] .= 0

    # Apply absorbing boundary conditions (Lysmer dampers)
    # For GPU: download small absorbing-node arrays, modify, scatter back
    if !(state.a isa Vector)
        a_abs_cpu = Array(state.a[absorbing_id])
        v_abs_cpu = Array(state.v[absorbing_id])
        for i in eachindex(absorbing_id)
            a_abs_cpu[i] -= absorb_matrix[i] * v_abs_cpu[i]
        end
        state.a[absorbing_id] .= CuArray(a_abs_cpu)
    else
        @inbounds @simd for i in eachindex(absorbing_id)
            state.a[absorbing_id[i]] -= absorb_matrix[i] * state.v[absorbing_id[i]]
        end
    end

    #--------------------------------------------------
    # Rate-state fault boundary
    #--------------------------------------------------

    # Step 4: Compute stick traction (free velocity if no friction)
    # FaultVFree = 2*v[n+1/2] + dt*a_internal/M
    # Download fault nodes from GPU for scalar computation
    v_fault_cpu  = Array(state.v[fault_id])
    a_fault_cpu  = Array(state.a[fault_id])
    M_fault_cpu  = Array(M_global[fault_id])
    fault_vfree_cpu = similar(M_fault_cpu)
    for i in eachindex(fault_id)
        fault_vfree_cpu[i] = 2 * v_fault_cpu[i] + dt * a_fault_cpu[i] / M_fault_cpu[i]
    end

    # Fault slip rate from previous timestep (for state update)
    v_prev_fault_cpu = Array(state.v_prev[fault_id])
    Vf_prev = 2 .* v_prev_fault_cpu .+ params.Vpl

    # Download fault state arrays for NR loop
    ψ_cpu  = Array(state.ψ)
    Vf_cpu = Array(state.Vf)
    τf_cpu = Array(state.τf)

    # Workspace for two-iteration NR search
    Vf_first_iter = similar(Vf_cpu)

    # Step 5-7: Solve nonlinear fault boundary equations (CPU loop)
    for i in eachindex(fault_id)
        # Save ψ before update so we can restore on NR failure
        ψ_saved = ψ_cpu[i]

        # Update state variable using slip rate from previous step
        ψ_cpu[i] = state_time_evolution(ψ_cpu[i], Vf_prev[i], dt,
                                        ics.friction.Lc[i], params.Vo)

        # Store slip rate for second iteration
        Vf_first_iter[i] = Vf_prev[i]

        # CRITICAL FIX: Convert perturbation stress → total stress for nr_search
        τ_total_guess = τf_cpu[i] + ics.τo[i]

        # Newton-Raphson search (1st iteration)
        Vf_cpu[i], τf_cpu[i] = nr_search(
            τ_total_guess,
            params.fo,
            params.Vo,
            ics.friction.a[i],
            ics.friction.b[i],
            ics.σo[i],
            ics.τo[i],
            ψ_cpu[i],
            fault_z[i],
            fault_vfree_cpu[i]
        )

        # Check for NaN (indicates NR failure) — restore ψ and continue
        if isnan(Vf_cpu[i]) || isnan(τf_cpu[i])
            ψ_cpu[i] = ψ_saved
            Vf_cpu[i] = Vf_prev[i]
            τf_cpu[i] = τ_total_guess - ics.τo[i]
            @warn "NR search produced NaN, restoring previous state" location=i iter=state.iteration
            continue
        end

        # 2nd iteration: Update state with average slip rate for accuracy
        ψ_after_first = ψ_cpu[i]
        avg_slip_rate = 0.5 * (Vf_cpu[i] + Vf_first_iter[i])
        ψ_cpu[i] = state_time_evolution(ψ_cpu[i], avg_slip_rate, dt,
                                        ics.friction.Lc[i], params.Vo)

        Vf_2nd, τf_2nd = nr_search(
            τf_cpu[i],
            params.fo,
            params.Vo,
            ics.friction.a[i],
            ics.friction.b[i],
            ics.σo[i],
            ics.τo[i],
            ψ_cpu[i],
            fault_z[i],
            fault_vfree_cpu[i]
        )

        if isnan(Vf_2nd) || isnan(τf_2nd)
            ψ_cpu[i] = ψ_after_first
            @warn "NR 2nd iteration produced NaN, keeping 1st iteration result" location=i
        else
            Vf_cpu[i] = Vf_2nd
            τf_cpu[i] = τf_2nd
        end
    end

    # Upload NR results back to state (works for both CPU and GPU)
    copyto!(state.ψ, _as_storage(state.ψ, ψ_cpu))
    copyto!(state.Vf, _as_storage(state.Vf, Vf_cpu))
    copyto!(state.fault_vfree, _as_storage(state.fault_vfree, fault_vfree_cpu))

    # Convert fault stress to perturbation from initial stress
    τf_pert_cpu = τf_cpu .- ics.τo
    copyto!(state.τf, _as_storage(state.τf, τf_pert_cpu))

    # Step 8: Apply fault traction to force vector
    # Download a[fault_id], modify, scatter back (for GPU path)
    if !(state.a isa Vector)
        a_flt_cpu = Array(state.a[fault_id])
        for i in eachindex(fault_id)
            a_flt_cpu[i] -= fault_matrix[i] * τf_pert_cpu[i]
        end
        state.a[fault_id] .= CuArray(a_flt_cpu)
    else
        @inbounds @simd for i in eachindex(fault_id)
            state.a[fault_id[i]] -= fault_matrix[i] * τf_pert_cpu[i]
        end
    end

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
