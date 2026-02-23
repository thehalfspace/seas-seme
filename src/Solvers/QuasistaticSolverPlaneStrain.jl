"""
    QuasistaticSolverPlaneStrain

Quasi-static earthquake cycle solver for plane-strain formulation.
AMG-preconditioned CG on 2*ndof system.

Same algorithm as antiplane (two-pass iteration, Kaneko et al. 2008/2011)
but operating on 2-component displacement/velocity fields.
"""

"""
    QuasistaticSolverPlaneStrain{P,O,T}

Quasi-static solver for plane-strain with AMG preconditioner.

# Fields
- `preconditioner::P`: AMG preconditioner (on 2*ndof reduced system)
- `stiffness_op::O`: StiffnessOperatorPlaneStrain
- `tolerance::T`: CG convergence tolerance
- `max_iterations::Int`: Maximum CG iterations
- `verbose::Bool`: Print convergence info
"""
struct QuasistaticSolverPlaneStrain{P,O,T<:AbstractFloat}
    preconditioner::P
    stiffness_op::O
    tolerance::T
    max_iterations::Int
    verbose::Bool
end


"""
    build_quasistatic_solver_plane_strain(weights, H, Ht, dof_id, mesh, fltni; ...)

Construct plane-strain quasi-static solver with AMG preconditioner.

# Arguments
- `weights::MetricWeightsPlaneStrain{T}`: Compact metric weight arrays
- `H::Matrix{T}`: Derivative matrix
- `Ht::Matrix{T}`: H' (pre-transposed)
- `dof_id`: DOF connectivity [nnodes, nnodes, n_elements] (spatial)
- `mesh`: UnstructuredSEMesh
- `fltni`: Free DOF indices in the 2*ndof space
- `tolerance`, `max_iterations`, `amg_max_levels`, `verbose`: Solver params

# Returns
- `QuasistaticSolverPlaneStrain` ready for time stepping
"""
function build_quasistatic_solver_plane_strain(
    weights::MetricWeightsPlaneStrain{T},
    H::Matrix{T},
    Ht::Matrix{T},
    dof_id::Array{Int,3},
    mesh,
    fltni::Vector{Int};
    tolerance::T=T(1e-6),
    max_iterations::Int=100,
    amg_max_levels::Int=10,
    verbose::Bool=true
) where T<:AbstractFloat

    ndof = mesh.ndof

    if verbose
        println("Building plane-strain AMG preconditioner (materializing K_el from metric weights)...")
    end

    # Assemble sparse stiffness (2*ndof × 2*ndof) — materializes K_el from weights
    K_sparse = stiffness_assembly_plane_strain(weights, H, Ht, dof_id, ndof)
    K_reduced = K_sparse[fltni, fltni]

    # Build AMG preconditioner
    ml = ruge_stuben(K_reduced, max_levels=amg_max_levels)
    amg_precond = aspreconditioner(ml)

    if verbose
        println("  AMG levels: ", length(ml.levels))
        println("  Coarsest level size: ", size(ml.levels[end].A, 1))
    end

    # Build matrix-free operator (tensor-product, zero-allocation per iteration)
    stiffness_op = StiffnessOperatorPlaneStrain(weights, H, Ht, dof_id,
                                                mesh.n_elements, ndof, fltni)

    return QuasistaticSolverPlaneStrain(amg_precond, stiffness_op, tolerance,
                                        max_iterations, verbose)
end


"""
    quasistatic_step!(state::SimulationStatePlaneStrain, solver, mesh, physics,
                      ics, params, dt)

Perform one quasi-static time step (plane-strain).

Same algorithm as antiplane (Kaneko et al. 2008/2011):
1. Store previous state
2. Two-pass iteration:
   a. Prescribe displacement on fault/creep boundaries
   b. Solve K*u = -f for interior DOFs
   c. Extract fault shear traction via tangential projection
   d. Update state variable and slip rate via rate-state friction
3. Update global velocity
"""
function quasistatic_step!(state::SimulationStatePlaneStrain{T},
                           solver::QuasistaticSolverPlaneStrain,
                           mesh, physics, ics, params, dt) where T

    fault_id = mesh.boundaries.fault.node_ids
    creep_id = mesh.boundaries.creep.node_ids
    fault_matrix = mesh.boundaries.fault.matrix
    mask = mesh.active_fault_mask
    ndof = state.ndof
    tangent = state.fault_tangent

    # Store previous state
    copyto!(state.u_prev, state.u)
    copyto!(state.v_prev, state.v)

    # Fault slip rate (total = 2*v_tang + Vpl)
    Vf_old = get_fault_tangential_velocity(state.v_prev, fault_id, tangent, ndof, params.Vpl)
    Vf_new = copy(Vf_old)

    if state.iteration == 0
        Vf_old .= get_fault_tangential_velocity(state.v, fault_id, tangent, ndof, params.Vpl)
        Vf_new .= copy(Vf_old)
    end

    op = solver.stiffness_op

    # Two-pass quasi-static iteration
    for pass in 1:2
        # Step 1: Prescribed displacement on boundaries
        fill!(state.f, zero(T))

        # Fault: u = u_prev + v * dt (both components)
        for i in eachindex(fault_id)
            nid = fault_id[i]
            state.f[nid]        = state.u_prev[nid]        + state.v[nid]        * dt
            state.f[ndof + nid] = state.u_prev[ndof + nid] + state.v[ndof + nid] * dt
        end

        # Creep: u = u_prev + v * dt (both components, v=0 for creep)
        for nid in creep_id
            state.f[nid]        = state.u_prev[nid]        + state.v[nid]        * dt
            state.f[ndof + nid] = state.u_prev[ndof + nid] + state.v[ndof + nid] * dt
        end

        # Step 2: Solve K*u = -f for free DOFs
        # Use operator's pre-allocated y_full buffer for K*f (rhs)
        apply_stiffness_plane_strain!(op.y_full, state.f,
                                      op.weights, op.H, op.Ht,
                                      op.dof_id, mesh.n_elements, ndof)
        rhs = op.y_full[op.fltni]

        # Warm-start CG from previous displacement (critical after dynamic→QS transition)
        x0 = state.u_prev[op.fltni]
        u_sol, history = cg!(x0, op, -rhs,
                            Pl=solver.preconditioner,
                            reltol=solver.tolerance,
                            maxiter=solver.max_iterations,
                            log=true)

        state.u[op.fltni] .= u_sol

        if any(isnan.(u_sol)) || any(isinf.(u_sol))
            @error "Non-finite values in plane-strain CG solution" pass iteration=state.iteration
            error("CG solver produced non-finite values")
        end

        # Enforce prescribed displacement on boundaries
        for i in eachindex(fault_id)
            nid = fault_id[i]
            state.u[nid]        = state.u_prev[nid]        + state.v[nid]        * dt
            state.u[ndof + nid] = state.u_prev[ndof + nid] + state.v[ndof + nid] * dt
        end
        for nid in creep_id
            state.u[nid]        = state.u_prev[nid]        + state.v[nid]        * dt
            state.u[ndof + nid] = state.u_prev[ndof + nid] + state.v[ndof + nid] * dt
        end

        # Step 3: Compute fault traction from K*u
        fill!(state.a, zero(T))
        apply_stiffness_plane_strain!(state.a, state.u,
                                      op.weights, op.H, op.Ht,
                                      op.dof_id, mesh.n_elements, ndof)

        # Zero force on creep boundary
        state.a[creep_id] .= 0
        state.a[ndof .+ creep_id] .= 0

        # Extract tangential shear traction and normal stress perturbation
        τf_new, σn_pert = fault_traction_from_force_plane_strain(
            state.a, fault_id, fault_matrix, tangent, state.fault_normal, ndof
        )
        state.τf .= τf_new
        state.σn_perturbation .= σn_pert

        # Steps 4 & 5: Update state variable and slip rate
        for i in eachindex(fault_id)
            if !mask[i]
                # Excluded endpoint: lock to plate rate, no state evolution
                Vf_new[i] = params.Vpl
                continue
            end

            state.ψ[i] = state_time_evolution(state.ψ[i], Vf_new[i], dt,
                                             ics.friction.Lc[i], params.Vo)

            # Normal stress: use constant σ_n⁰ for now.
            # TODO: For dipping faults, use σn_eff = ics.σo[i] + state.σn_perturbation[i]
            # The perturbation extraction needs calibration (subtract initial state contribution)
            # before it can be used in the friction law.
            Vf_new[i] = fault_slip_rate(state.ψ[i], state.τf[i],
                                       ics.τo[i], ics.σo[i],
                                       ics.friction.a[i], ics.friction.b[i],
                                       params.Vo, params.fo)

            if !isfinite(Vf_new[i]) && state.iteration == 1
                @error "Non-finite slip rate at fault node" i pass
                error("Non-finite slip rate computed")
            end
        end

        # Average slip rates
        Vf_new .= 0.5 .* (Vf_old .+ Vf_new)

        # Update fault velocity (projected onto tangent)
        set_fault_velocity_plane_strain!(state.v, fault_id, Vf_new, params.Vpl,
                                         tangent, ndof)

        # Zero creep velocity
        state.v[creep_id] .= 0
        state.v[ndof .+ creep_id] .= 0
    end  # Two-pass iteration

    # Update global velocity from displacement change
    state.v .= (state.u .- state.u_prev) ./ dt

    # Safety: clamp velocity to physical bounds (shear wave velocity)
    # Prevents unphysical spikes from CG errors amplified by small dt
    clamp!(state.v, -physics.vs, physics.vs)

    # Re-enforce boundary velocities
    set_fault_velocity_plane_strain!(state.v, fault_id, Vf_new, params.Vpl,
                                     tangent, ndof)
    state.v[creep_id] .= 0
    state.v[ndof .+ creep_id] .= 0

    # Blowup detection
    v_max = maximum(abs.(state.v))
    u_max = maximum(abs.(state.u))
    if v_max > 1e10 || u_max > 1e10
        @error "PS-QS blowup detected" v_max u_max dt iteration=state.iteration Vf_max=maximum(abs.(Vf_new))
        error("Quasistatic solver produced non-physical values")
    end

    # Zero prescribed boundary accelerations
    fill!(state.a, zero(T))
    state.u[creep_id] .= 0
    state.u[ndof .+ creep_id] .= 0

    # Store fault slip rate
    state.Vf .= Vf_new

    return nothing
end


# ============================================================================
# GPU quasi-static solver (plane-strain) using AMGX PCG+AMG
# ============================================================================

"""
    build_quasistatic_solver_gpu_plane_strain(weights, H, Ht, dof_id, mesh, fltni; ...)

Build a GPU quasi-static solver for plane-strain using AMGX PCG+AMG.
"""
function build_quasistatic_solver_gpu_plane_strain(
    weights::MetricWeightsPlaneStrain{T},
    H::Matrix{T},
    Ht::Matrix{T},
    dof_id::Array{Int,3},
    mesh,
    fltni::Vector{Int};
    tolerance::T=T(1e-8),
    amgx_config_str::String=DEFAULT_AMGX_CONFIG
) where T<:AbstractFloat

    println("Building plane-strain GPU AMGX PCG+AMG solver...")
    ndof = mesh.ndof

    K_sparse = stiffness_assembly_plane_strain(weights, H, Ht, dof_id, ndof)
    K_reduced = K_sparse[fltni, fltni]

    return _build_amgx_solver(K_reduced, fltni, tolerance, amgx_config_str)
end


"""
    quasistatic_step!(state::SimulationStatePlaneStrain, solver::QuasistaticSolverGPU, ...)

GPU quasi-static step for plane-strain using AMGX.
Fault NR loop runs on CPU; matvec runs on GPU via CuArray dispatch.
"""
function quasistatic_step!(state::SimulationStatePlaneStrain{T},
                           solver::QuasistaticSolverGPU{T},
                           mesh, physics, ics, params, dt,
                           weights_gpu, H_gpu, Ht_gpu, dof_id_gpu, M_global_gpu) where T

    fault_id = mesh.boundaries.fault.node_ids
    creep_id = mesh.boundaries.creep.node_ids
    fault_matrix = mesh.boundaries.fault.matrix
    mask = mesh.active_fault_mask
    ndof = state.ndof
    tangent = state.fault_tangent   # CPU matrix [2, nfault]

    copyto!(state.u_prev, state.u)
    copyto!(state.v_prev, state.v)

    # Compute Vf_old on CPU (tangential velocity projection, fault nodes only)
    v_prev_cpu = Array(state.v_prev)
    Vf_old = get_fault_tangential_velocity(v_prev_cpu, fault_id, tangent, ndof, params.Vpl)
    Vf_new = copy(Vf_old)

    if state.iteration == 0
        v_cpu0 = Array(state.v)
        Vf_old .= get_fault_tangential_velocity(v_cpu0, fault_id, tangent, ndof, params.Vpl)
        Vf_new .= copy(Vf_old)
    end

    fltni_gpu = solver.fltni_gpu

    # Pre-download boundary arrays once (small, used in every pass)
    u_prev_fault_x = Array(state.u_prev[fault_id])
    u_prev_fault_y = Array(state.u_prev[ndof .+ fault_id])
    u_prev_creep_x = Array(state.u_prev[creep_id])
    u_prev_creep_y = Array(state.u_prev[ndof .+ creep_id])

    for pass in 1:2
        # Step 1: Prescribed displacement on boundaries.
        # Download v at boundary nodes, compute prescribed u, scatter-upload.
        fill!(state.f, zero(T))
        v_fault_x = Array(state.v[fault_id])
        v_fault_y = Array(state.v[ndof .+ fault_id])
        v_creep_x = Array(state.v[creep_id])
        v_creep_y = Array(state.v[ndof .+ creep_id])

        f_fault_x = u_prev_fault_x .+ v_fault_x .* dt
        f_fault_y = u_prev_fault_y .+ v_fault_y .* dt
        f_creep_x = u_prev_creep_x .+ v_creep_x .* dt
        f_creep_y = u_prev_creep_y .+ v_creep_y .* dt

        state.f[fault_id]        .= CuArray(f_fault_x)
        state.f[ndof .+ fault_id] .= CuArray(f_fault_y)
        state.f[creep_id]        .= CuArray(f_creep_x)
        state.f[ndof .+ creep_id] .= CuArray(f_creep_y)

        # Step 2: K*f (GPU kernel) → RHS for AMGX
        apply_stiffness_plane_strain!(state.a, state.f, weights_gpu, H_gpu, Ht_gpu,
                                      dof_id_gpu, mesh.n_elements, ndof)
        rhs_gpu = Array(state.a[fltni_gpu])

        x0_gpu = Array(state.u_prev[fltni_gpu])
        AMGX.upload!(solver.b_buf, -rhs_gpu)
        AMGX.upload!(solver.x_buf, x0_gpu)
        AMGX.solve!(solver.x_buf, solver.amgx_solver, solver.b_buf)

        u_sol_cpu = AMGX.download(solver.x_buf)
        if any(isnan.(u_sol_cpu)) || any(isinf.(u_sol_cpu))
            @error "Non-finite values in PS AMGX solution" pass iteration=state.iteration
            error("AMGX solver produced non-finite values")
        end
        state.u[fltni_gpu] .= CuArray(u_sol_cpu)

        # Enforce prescribed displacements (scatter-upload)
        state.u[fault_id]        .= CuArray(f_fault_x)
        state.u[ndof .+ fault_id] .= CuArray(f_fault_y)
        state.u[creep_id]        .= CuArray(f_creep_x)
        state.u[ndof .+ creep_id] .= CuArray(f_creep_y)

        # Step 3: Compute fault traction from K*u (GPU kernel)
        fill!(state.a, zero(T))
        apply_stiffness_plane_strain!(state.a, state.u, weights_gpu, H_gpu, Ht_gpu,
                                      dof_id_gpu, mesh.n_elements, ndof)
        state.a[creep_id] .= 0
        state.a[ndof .+ creep_id] .= 0

        # Extract traction on CPU
        a_cpu = Array(state.a)
        τf_new, σn_pert = fault_traction_from_force_plane_strain(
            a_cpu, fault_id, fault_matrix, tangent, state.fault_normal, ndof
        )
        copyto!(state.τf, CuArray(τf_new))
        copyto!(state.σn_perturbation, CuArray(σn_pert))

        # Steps 4 & 5: Rate-state evolution (CPU)
        ψ_cpu = Array(state.ψ)

        for i in eachindex(fault_id)
            if !mask[i]
                Vf_new[i] = params.Vpl
                continue
            end

            ψ_cpu[i] = state_time_evolution(ψ_cpu[i], Vf_new[i], dt,
                                            ics.friction.Lc[i], params.Vo)
            Vf_new[i] = fault_slip_rate(ψ_cpu[i], τf_new[i],
                                        ics.τo[i], ics.σo[i],
                                        ics.friction.a[i], ics.friction.b[i],
                                        params.Vo, params.fo)

            if !isfinite(Vf_new[i]) && state.iteration == 1
                @error "Non-finite slip rate at fault node" i pass
                error("Non-finite slip rate computed")
            end
        end

        copyto!(state.ψ, CuArray(ψ_cpu))

        # Average slip rates
        Vf_new .= 0.5 .* (Vf_old .+ Vf_new)

        # Update fault velocity: project Vf onto tangent, upload full v vector
        v_cpu_full = Array(state.v)
        set_fault_velocity_plane_strain!(v_cpu_full, fault_id, Vf_new, params.Vpl, tangent, ndof)
        for nid in creep_id
            v_cpu_full[nid]        = zero(T)
            v_cpu_full[ndof + nid] = zero(T)
        end
        copyto!(state.v, CuArray(v_cpu_full))
    end  # Two-pass iteration

    # Update velocity from displacement change
    state.v .= (state.u .- state.u_prev) ./ dt
    clamp!(state.v, -physics.vs, physics.vs)

    # Re-enforce boundary velocities (download full v, modify, upload)
    v_cpu_final = Array(state.v)
    set_fault_velocity_plane_strain!(v_cpu_final, fault_id, Vf_new, params.Vpl, tangent, ndof)
    for nid in creep_id
        v_cpu_final[nid]        = zero(T)
        v_cpu_final[ndof + nid] = zero(T)
    end
    copyto!(state.v, CuArray(v_cpu_final))

    v_max = maximum(abs.(Array(state.v)))
    u_max = maximum(abs.(Array(state.u)))
    if v_max > 1e10 || u_max > 1e10
        @error "PS-GPU-QS blowup detected" v_max u_max dt iteration=state.iteration
        error("GPU plane-strain quasistatic solver produced non-physical values")
    end

    fill!(state.a, zero(T))
    state.u[creep_id]        .= zero(T)
    state.u[ndof .+ creep_id] .= zero(T)
    state.v[creep_id]        .= zero(T)
    state.v[ndof .+ creep_id] .= zero(T)

    copyto!(state.Vf, CuArray(Vf_new))

    return nothing
end
