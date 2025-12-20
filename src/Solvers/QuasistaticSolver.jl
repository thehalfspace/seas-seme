"""
    QuasistaticSolver

Quasi-static earthquake cycle solver with AMG-preconditioned conjugate gradient.

Solves the equilibrium equation K*u = -f iteratively using matrix-free operator
and algebraic multigrid preconditioner for efficient solution of large systems.

Based on Kaneko et al. (2008, 2011) formulation.
"""

# Note: All dependencies (IterativeSolvers, AlgebraicMultigrid, etc.) are imported
# at the main module level (SEAS_SEME.jl). Functions from other submodules
# (Discretization, Physics) are available in the same namespace.


"""
    QuasistaticSolver{P,O,T}

Quasi-static solver with AMG preconditioner and matrix-free operator.

# Fields
- `preconditioner::P`: AMG preconditioner (AlgebraicMultigrid.MultiLevel)
- `stiffness_op::O`: Matrix-free stiffness operator
- `tolerance::T`: Relative tolerance for CG convergence
- `max_iterations::Int`: Maximum CG iterations
- `verbose::Bool`: Print convergence information
"""
struct QuasistaticSolver{P,O,T<:AbstractFloat}
    preconditioner::P
    stiffness_op::O
    tolerance::T
    max_iterations::Int
    verbose::Bool
end


"""
    build_quasistatic_solver(K_el, dof_id, mesh, fltni;
                            tolerance=1e-6, max_iterations=100,
                            amg_max_levels=10, verbose=true)

Construct quasi-static solver with AMG preconditioner.

# Arguments
- `K_el`: Elemental stiffness matrices [nnodes², nnodes², n_elements]
- `dof_id`: DOF connectivity [nnodes, nnodes, n_elements]
- `mesh`: UnstructuredSEMesh object
- `fltni`: Non-fault/creep DOF indices (free DOFs)
- `tolerance`: CG convergence tolerance (default: 1e-6)
- `max_iterations`: Maximum CG iterations (default: 100)
- `amg_max_levels`: Maximum levels in AMG hierarchy (default: 10)
- `verbose`: Print convergence info (default: true)

# Returns
- `QuasistaticSolver` instance ready for time stepping

# Notes
- Assembles sparse stiffness matrix once for AMG preconditioner setup
- Creates matrix-free operator for efficient CG iterations
- AMG preconditioner significantly accelerates convergence for high-order elements
"""
function build_quasistatic_solver(K_el::Array{T,3}, dof_id::Array{Int,3},
                                  mesh, fltni::Vector{Int};
                                  tolerance::T=T(1e-6),
                                  max_iterations::Int=100,
                                  amg_max_levels::Int=10,
                                  verbose::Bool=true) where T<:AbstractFloat

    # Build sparse stiffness matrix for preconditioner
    # (only need to do this once - amortized cost)
    if verbose
        println("Building AMG preconditioner...")
    end

    K_sparse = stiffness_assembly(K_el, dof_id)
    K_reduced = K_sparse[fltni, fltni]

    # Build AMG preconditioner (Ruge-Stuben)
    ml = ruge_stuben(K_reduced, max_levels=amg_max_levels)
    amg_precond = aspreconditioner(ml)  # Wrap for IterativeSolvers.jl

    if verbose
        println("  AMG levels: ", length(ml.levels))
        println("  Coarsest level size: ", size(ml.levels[end].A, 1))
    end

    # Build matrix-free operator
    stiffness_op = StiffnessOperator(K_el, dof_id, mesh.n_elements,
                                    mesh.ndof, fltni)

    return QuasistaticSolver(amg_precond, stiffness_op, tolerance,
                            max_iterations, verbose)
end


"""
    quasistatic_step!(state, solver, mesh, physics, ics, params, dt)

Perform one quasi-static time step.

# Algorithm (from Kaneko et al. 2008, 2011)
1. Store previous displacement and velocity
2. Predict fault slip rate from previous state
3. Two-pass iteration:
   a. Compute prescribed displacement on fault/creep boundaries
   b. Solve K*u = -f for interior DOFs (CG with AMG)
   c. Update fault stress from K*u
   d. Update state variable ψ
   e. Compute new slip rate from rate-state friction
   f. Average slip rates between passes
4. Update global velocity from displacement change

# Arguments
- `state`: SimulationState (modified in place)
- `solver`: QuasistaticSolver
- `mesh`: UnstructuredSEMesh
- `physics`: Material properties
- `ics`: Initial conditions (fault parameters)
- `params`: Simulation parameters (Vpl, V₀, f₀)
- `dt`: Time step [s]

# Modifies
- `state.u`: Displacement field
- `state.v`: Velocity field
- `state.τf`: Fault shear stress
- `state.ψ`: Fault state variable
- `state.Vf`: Fault slip rate
"""
function quasistatic_step!(state, solver::QuasistaticSolver, mesh, physics,
                          ics, params, dt)
    # Extract boundary indices
    fault_id = mesh.boundaries.fault.node_ids
    creep_id = mesh.boundaries.creep.node_ids
    fault_matrix = mesh.boundaries.fault.matrix

    # Store previous state
    copyto!(state.u_prev, state.u)
    copyto!(state.v_prev, state.v)

    # Fault slip rate (total velocity = 2*v + Vpl for single-sided fault)
    Vf_old = 2 * state.v_prev[fault_id] .+ params.Vpl
    Vf_new = copy(Vf_old)

    # Special handling for first iteration: use current state.v instead of v_prev
    if state.iteration == 0
        Vf_old .= 2 * state.v[fault_id] .+ params.Vpl
        Vf_new .= copy(Vf_old)
    end

    # Two-pass quasi-static iteration for accuracy
    for pass in 1:2

        # Step 1: Prescribed displacement on boundaries
        fill!(state.f, zero(eltype(state.f)))
        state.f[fault_id] .= state.u_prev[fault_id] .+ state.v[fault_id] .* dt
        state.f[creep_id] .= state.u_prev[creep_id] .+ state.v[creep_id] .* dt

        # Step 2: Solve K*u = -f for free DOFs (CG with AMG preconditioner)
        rhs_full = apply_stiffness!(zeros(eltype(state.u), mesh.ndof),
                                   state.f, solver.stiffness_op.K_el,
                                   solver.stiffness_op.dof_id,
                                   mesh.n_elements)
        rhs = rhs_full[solver.stiffness_op.fltni]

        # Conjugate gradient with AMG preconditioner
        u_sol, history = cg(solver.stiffness_op, -rhs,
                           Pl=solver.preconditioner,
                           reltol=solver.tolerance,
                           maxiter=solver.max_iterations,
                           log=true)

        state.u[solver.stiffness_op.fltni] .= u_sol

        # Convergence diagnostics
        if solver.verbose && (state.iteration <= 10 || mod(state.iteration, 500) == 0)
            niter = length(history[:resnorm]) - 1
            converged = history.isconverged
            println("  QS CG iterations: $niter, converged: $converged")
            if !converged
                println("    Final residual: $(history[:resnorm][end])")
                println("    u_sol stats: min=$(minimum(u_sol)), max=$(maximum(u_sol)), any_nan=$(any(isnan.(u_sol)))")
            end
        end

        # Check for NaN in solution
        if any(isnan.(u_sol)) || any(isinf.(u_sol))
            @error "Non-finite values in CG solution" pass iteration=state.iteration converged=history.isconverged
            println("  rhs stats: min=$(minimum(rhs)), max=$(maximum(rhs)), any_nan=$(any(isnan.(rhs)))")
            println("  u_sol stats: min=$(minimum(u_sol)), max=$(maximum(u_sol)), any_nan=$(any(isnan.(u_sol)))")
            error("CG solver produced non-finite values")
        end

        # Enforce prescribed displacement on boundaries
        state.u[fault_id] .= state.u_prev[fault_id] .+ state.v[fault_id] .* dt
        state.u[creep_id] .= state.u_prev[creep_id] .+ state.v[creep_id] .* dt

        # Step 3: Compute fault traction
        fill!(state.a, zero(eltype(state.a)))
        state.a .= apply_stiffness!(state.a, state.u,
                                   solver.stiffness_op.K_el,
                                   solver.stiffness_op.dof_id,
                                   mesh.n_elements)

        # Enforce zero traction on creep boundary
        state.a[creep_id] .= 0

        # Fault shear stress (τ = -K*u / fault_matrix)
        # Check for zeros in fault_matrix which would cause NaN
        if state.iteration == 1 && pass == 1
            n_zeros = sum(fault_matrix .== 0)
            if n_zeros > 0
                @error "Zero values in fault impedance matrix" n_zeros total=length(fault_matrix)
                println("  fault_matrix stats: min=$(minimum(fault_matrix)), max=$(maximum(fault_matrix))")
                println("  Zero indices: ", findall(fault_matrix .== 0))
                error("Cannot compute fault stress with zero impedance matrix values")
            end
        end

        state.τf .= -state.a[fault_id] ./ fault_matrix

        # Check for NaN in fault stress
        if any(isnan.(state.τf)) || any(isinf.(state.τf))
            @error "Non-finite fault stress" pass iteration=state.iteration
            println("  state.a[fault_id] stats: min=$(minimum(state.a[fault_id])), max=$(maximum(state.a[fault_id]))")
            println("  fault_matrix stats: min=$(minimum(fault_matrix)), max=$(maximum(fault_matrix))")
            println("  τf stats: min=$(minimum(state.τf)), max=$(maximum(state.τf)), any_nan=$(any(isnan.(state.τf)))")
            error("Non-finite fault stress computed")
        end

        # Steps 4 & 5: Update state variable and compute slip rate
        for i in eachindex(fault_id)
            # Update state variable (aging law)
            state.ψ[i] = state_time_evolution(state.ψ[i], Vf_new[i], dt,
                                             ics.friction.Lc[i], params.Vo)

            # Compute slip rate from rate-state friction
            Vf_new[i] = fault_slip_rate(state.ψ[i], state.τf[i],
                                       ics.τo[i], ics.σo[i],
                                       ics.friction.a[i], ics.friction.b[i],
                                       params.Vo, params.fo)

            # Check for NaN in slip rate
            if !isfinite(Vf_new[i]) && state.iteration == 1
                @error "Non-finite slip rate at fault node" i pass
                println("  ψ=$(state.ψ[i]), τf=$(state.τf[i]), τo=$(ics.τo[i]), σo=$(ics.σo[i])")
                println("  a=$(ics.friction.a[i]), b=$(ics.friction.b[i]), Vo=$(params.Vo), fo=$(params.fo)")
                error("Non-finite slip rate computed")
            end
        end

        # Average slip rates between old and new predictions
        Vf_new .= 0.5 .* (Vf_old .+ Vf_new)

        # Update fault velocity (single-sided fault: v = 0.5*(Vf - Vpl))
        state.v[fault_id] .= 0.5 .* (Vf_new .- params.Vpl)
        state.v[creep_id] .= 0
    end  # Two-pass iteration

    # Update global velocity from displacement change
    state.v .= (state.u .- state.u_prev) ./ dt
    state.v[fault_id] .= 0.5 .* (Vf_new .- params.Vpl)
    state.v[creep_id] .= 0

    # Zero out prescribed boundary accelerations
    fill!(state.a, zero(eltype(state.a)))
    state.u[creep_id] .= 0
    state.v[creep_id] .= 0

    # Store fault slip rate for output
    state.Vf .= Vf_new

    return nothing
end
