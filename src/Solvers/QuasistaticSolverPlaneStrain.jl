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
    build_quasistatic_solver_plane_strain(K_el, dof_id, mesh, fltni;
                                          tolerance=1e-6, max_iterations=100,
                                          amg_max_levels=10, verbose=true)

Construct plane-strain quasi-static solver with AMG preconditioner.

# Arguments
- `K_el`: Elemental stiffness matrices [2*N, 2*N, n_elements]
- `dof_id`: DOF connectivity [nnodes, nnodes, n_elements] (spatial)
- `mesh`: UnstructuredSEMesh
- `fltni`: Free DOF indices in the 2*ndof space
- `tolerance`, `max_iterations`, `amg_max_levels`, `verbose`: Solver params

# Returns
- `QuasistaticSolverPlaneStrain` ready for time stepping
"""
function build_quasistatic_solver_plane_strain(
    K_el::Array{T,3}, dof_id::Array{Int,3},
    mesh, fltni::Vector{Int};
    tolerance::T=T(1e-6),
    max_iterations::Int=100,
    amg_max_levels::Int=10,
    verbose::Bool=true
) where T<:AbstractFloat

    ndof = mesh.ndof

    if verbose
        println("Building plane-strain AMG preconditioner...")
    end

    # Assemble sparse stiffness (2*ndof × 2*ndof)
    K_sparse = stiffness_assembly_plane_strain(K_el, dof_id, ndof)
    K_reduced = K_sparse[fltni, fltni]

    # Build AMG preconditioner
    ml = ruge_stuben(K_reduced, max_levels=amg_max_levels)
    amg_precond = aspreconditioner(ml)

    if verbose
        println("  AMG levels: ", length(ml.levels))
        println("  Coarsest level size: ", size(ml.levels[end].A, 1))
        println("  Total DOFs (2*ndof): ", 2 * ndof)
        println("  Free DOFs: ", length(fltni))
    end

    # Build matrix-free operator
    stiffness_op = StiffnessOperatorPlaneStrain(K_el, dof_id, mesh.n_elements,
                                                 ndof, fltni)

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
        rhs_full = zeros(T, 2 * ndof)
        apply_stiffness_plane_strain!(rhs_full, state.f, solver.stiffness_op.K_el,
                                      solver.stiffness_op.dof_id, mesh.n_elements, ndof)
        rhs = rhs_full[solver.stiffness_op.fltni]

        u_sol, history = cg(solver.stiffness_op, -rhs,
                           Pl=solver.preconditioner,
                           reltol=solver.tolerance,
                           maxiter=solver.max_iterations,
                           log=true)

        state.u[solver.stiffness_op.fltni] .= u_sol

        if solver.verbose && (state.iteration <= 10 || mod(state.iteration, 500) == 0)
            niter = length(history[:resnorm]) - 1
            converged = history.isconverged
            println("  PS-QS CG iterations: $niter, converged: $converged")
        end

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
        apply_stiffness_plane_strain!(state.a, state.u, solver.stiffness_op.K_el,
                                      solver.stiffness_op.dof_id, mesh.n_elements, ndof)

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

    # Re-enforce boundary velocities
    set_fault_velocity_plane_strain!(state.v, fault_id, Vf_new, params.Vpl,
                                     tangent, ndof)
    state.v[creep_id] .= 0
    state.v[ndof .+ creep_id] .= 0

    # Zero prescribed boundary accelerations
    fill!(state.a, zero(T))
    state.u[creep_id] .= 0
    state.u[ndof .+ creep_id] .= 0

    # Store fault slip rate
    state.Vf .= Vf_new

    return nothing
end
