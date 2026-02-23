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
    build_quasistatic_solver(weights, H, Ht, dof_id, mesh, fltni; ...)

Construct quasi-static solver with AMG preconditioner.

# Arguments
- `weights::MetricWeightsAntiplane{T}`: Compact metric weight arrays
- `H::Matrix{T}`: Derivative matrix
- `Ht::Matrix{T}`: H' (pre-transposed)
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
- Assembles sparse stiffness matrix once (materializes K_el from weights) for AMG setup
- Creates matrix-free operator using tensor-product matvec for efficient CG iterations
- AMG preconditioner significantly accelerates convergence for high-order elements
"""
function build_quasistatic_solver(
    weights::MetricWeightsAntiplane{T},
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

    # Build sparse stiffness matrix for preconditioner
    # (only need to do this once - amortized cost; materializes K_el from weights)
    if verbose
        println("Building AMG preconditioner (materializing K_el from metric weights)...")
    end

    K_sparse = stiffness_assembly(weights, H, Ht, dof_id)
    K_reduced = K_sparse[fltni, fltni]

    # Build AMG preconditioner (Ruge-Stuben)
    ml = ruge_stuben(K_reduced, max_levels=amg_max_levels)
    amg_precond = aspreconditioner(ml)  # Wrap for IterativeSolvers.jl

    if verbose
        println("  AMG levels: ", length(ml.levels))
        println("  Coarsest level size: ", size(ml.levels[end].A, 1))
    end

    # Build matrix-free operator (tensor-product, zero-allocation per iteration)
    stiffness_op = StiffnessOperator(weights, H, Ht, dof_id, mesh.n_elements,
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
        # Use operator's y_full workspace to compute K*f (rhs)
        op = solver.stiffness_op
        apply_stiffness!(op.y_full, state.f, op.weights, op.H, op.Ht, op.dof_id, mesh.n_elements)
        rhs = op.y_full[op.fltni]

        # Warm-start CG from previous displacement (critical after dynamic→QS transition)
        x0 = state.u_prev[op.fltni]
        u_sol, history = cg!(x0, op, -rhs,
                            Pl=solver.preconditioner,
                            reltol=solver.tolerance,
                            maxiter=solver.max_iterations,
                            log=true)

        state.u[op.fltni] .= u_sol

        # Check for NaN in solution
        if any(isnan.(u_sol)) || any(isinf.(u_sol))
            @error "Non-finite values in CG solution" pass iteration=state.iteration converged=history.isconverged
            error("CG solver produced non-finite values")
        end

        # Enforce prescribed displacement on boundaries
        state.u[fault_id] .= state.u_prev[fault_id] .+ state.v[fault_id] .* dt
        state.u[creep_id] .= state.u_prev[creep_id] .+ state.v[creep_id] .* dt

        # Step 3: Compute fault traction
        fill!(state.a, zero(eltype(state.a)))
        apply_stiffness!(state.a, state.u, op.weights, op.H, op.Ht, op.dof_id, mesh.n_elements)

        # Enforce zero traction on creep boundary
        state.a[creep_id] .= 0

        # Fault shear stress (τ = -K*u / fault_matrix)
        state.τf .= -state.a[fault_id] ./ fault_matrix

        if any(isnan.(state.τf)) || any(isinf.(state.τf))
            @error "Non-finite fault stress" pass iteration=state.iteration
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

    # Safety: clamp velocity to physical bounds (shear wave velocity)
    # Prevents unphysical spikes from CG errors amplified by small dt
    clamp!(state.v, -physics.vs, physics.vs)

    # Re-enforce boundary velocities
    state.v[fault_id] .= 0.5 .* (Vf_new .- params.Vpl)
    state.v[creep_id] .= 0

    # Blowup detection
    v_max = maximum(abs.(state.v))
    u_max = maximum(abs.(state.u))
    if v_max > 1e10 || u_max > 1e10
        @error "QS blowup detected" v_max u_max dt iteration=state.iteration
        error("Quasistatic solver produced non-physical values")
    end

    # Zero out prescribed boundary accelerations
    fill!(state.a, zero(eltype(state.a)))
    state.u[creep_id] .= 0
    state.v[creep_id] .= 0

    # Store fault slip rate for output
    state.Vf .= Vf_new

    return nothing
end


# ============================================================================
# GPU quasi-static solver using AMGX PCG+AMG
# ============================================================================

"""
    QuasistaticSolverGPU{T}

GPU quasi-static solver using AMGX PCG+AMG, all on GPU.

# Fields
- `amgx_solver`: AMGX solver object (PCG with AMG preconditioner)
- `amgx_matrix`: AMGX matrix wrapping the sparse stiffness
- `amgx_resources`: AMGX resources (GPU memory pool)
- `amgx_config`: AMGX configuration
- `x_buf`: AMGX vector for solution
- `b_buf`: AMGX vector for RHS
- `fltni_gpu`: Free DOF indices on GPU
- `ndof_reduced`: Number of free DOFs
- `tolerance::T`: Solver tolerance
"""
struct QuasistaticSolverGPU{T<:AbstractFloat}
    amgx_solver
    amgx_matrix
    amgx_resources
    amgx_config
    x_buf
    b_buf
    fltni_gpu::CuVector{Int32}
    ndof_reduced::Int
    tolerance::T
end


"""
    build_quasistatic_solver_gpu(weights, H, Ht, dof_id, mesh, fltni; tolerance, amgx_config_str)

Build a GPU quasi-static solver using AMGX PCG+AMG on GPU.

Assembles the sparse stiffness matrix on CPU, uploads to GPU, and sets up AMGX.
"""
function build_quasistatic_solver_gpu(
    weights::MetricWeightsAntiplane{T},
    H::Matrix{T},
    Ht::Matrix{T},
    dof_id::Array{Int,3},
    mesh,
    fltni::Vector{Int};
    tolerance::T=T(1e-8),
    amgx_config_str::String=DEFAULT_AMGX_CONFIG
) where T<:AbstractFloat

    println("Building GPU AMGX PCG+AMG solver...")

    # Assemble sparse stiffness on CPU (one-time cost)
    K_sparse = stiffness_assembly(weights, H, Ht, dof_id)
    K_reduced = K_sparse[fltni, fltni]

    return _build_amgx_solver(K_reduced, fltni, tolerance, amgx_config_str)
end


const DEFAULT_AMGX_CONFIG = """{
    "config_version": 2,
    "solver": {
        "solver": "PCG",
        "preconditioner": {
            "solver": "AMG",
            "algorithm": "AGGREGATION",
            "cycle": "V",
            "smoother": "BLOCK_JACOBI",
            "presweeps": 1,
            "postsweeps": 1,
            "max_levels": 10,
            "coarse_solver": "DENSE_LU_SOLVER"
        },
        "max_iters": 200,
        "tolerance": 1e-8,
        "monitor_residual": 1,
        "convergence": "RELATIVE_INI_CORE"
    }
}"""


"""
    _build_amgx_solver(K_reduced, fltni, tolerance, amgx_config_str)

Internal: initialize AMGX resources, upload matrix, build solver.
"""
function _build_amgx_solver(K_reduced, fltni::Vector{Int}, tolerance::T, amgx_config_str::String) where T
    ndof_reduced = size(K_reduced, 1)

    # Initialize AMGX
    amgx_config    = AMGX.Config(amgx_config_str)
    amgx_resources = AMGX.Resources(amgx_config)

    # Upload sparse matrix to GPU (AMGX.jl handles 0-indexing internally)
    amgx_matrix = AMGX.AMGXMatrix(amgx_resources, AMGX.dDDI)
    K_cu = CUDA.CUSPARSE.CuSparseMatrixCSR(K_reduced)
    AMGX.upload!(amgx_matrix, K_cu)

    # Create and setup solver (AMG hierarchy built here)
    amgx_solver = AMGX.Solver(amgx_resources, AMGX.dDDI, amgx_config)
    AMGX.setup!(amgx_solver, amgx_matrix)

    # Allocate solution and RHS vectors on GPU
    x_buf = AMGX.AMGXVector(amgx_resources, AMGX.dDDI)
    b_buf = AMGX.AMGXVector(amgx_resources, AMGX.dDDI)
    AMGX.upload!(x_buf, zeros(T, ndof_reduced))
    AMGX.upload!(b_buf, zeros(T, ndof_reduced))

    fltni_gpu = CuArray(Int32.(fltni))

    println("  AMGX solver initialized ($(ndof_reduced) free DOFs)")

    return QuasistaticSolverGPU{T}(
        amgx_solver, amgx_matrix, amgx_resources, amgx_config,
        x_buf, b_buf, fltni_gpu, ndof_reduced, tolerance
    )
end


"""
    quasistatic_step!(state, solver::QuasistaticSolverGPU, mesh, physics, ics, params, dt,
                      weights_gpu, H_gpu, Ht_gpu, dof_id_gpu, M_global_gpu)

GPU quasi-static step using AMGX PCG+AMG.

Key differences from CPU path:
- CG replaced by AMGX.solve! (PCG+AMG fully on GPU)
- Fault NR loop runs on CPU: download ~1-2k fault arrays, run NR, upload back
- apply_stiffness! dispatches to GPU kernel via CuArray type
"""
function quasistatic_step!(state, solver::QuasistaticSolverGPU{T}, mesh, physics,
                           ics, params, dt,
                           weights_gpu, H_gpu, Ht_gpu, dof_id_gpu, M_global_gpu) where T
    fault_id = mesh.boundaries.fault.node_ids
    creep_id = mesh.boundaries.creep.node_ids
    fault_matrix = mesh.boundaries.fault.matrix
    ndof = mesh.ndof

    # Store previous state (GPU broadcast, zero-alloc)
    copyto!(state.u_prev, state.u)
    copyto!(state.v_prev, state.v)

    # Compute Vf_old on CPU (download fault-only arrays)
    v_fault_cpu = Array(state.v_prev[fault_id])
    Vf_old = 2 .* v_fault_cpu .+ params.Vpl
    Vf_new = copy(Vf_old)

    if state.iteration == 0
        v_fault_cpu0 = Array(state.v[fault_id])
        Vf_old .= 2 .* v_fault_cpu0 .+ params.Vpl
        Vf_new .= copy(Vf_old)
    end

    fltni_gpu = solver.fltni_gpu

    # Pre-download small boundary arrays to CPU (used in every pass)
    u_prev_fault_cpu = Array(state.u_prev[fault_id])
    u_prev_creep_cpu = Array(state.u_prev[creep_id])
    v_fault_cpu_init = Array(state.v[fault_id])
    v_creep_cpu_init = Array(state.v[creep_id])

    for pass in 1:2
        # Step 1: Prescribed displacement on boundaries.
        # Download v at boundary nodes, compute prescribed u, upload back.
        fill!(state.f, zero(T))
        v_fault_now = Array(state.v[fault_id])
        v_creep_now = Array(state.v[creep_id])

        # Compute prescribed displacements on CPU, scatter-upload to GPU
        f_fault_cpu = u_prev_fault_cpu .+ v_fault_now .* dt
        f_creep_cpu = u_prev_creep_cpu .+ v_creep_now .* dt
        state.f[fault_id] .= CuArray(f_fault_cpu)
        state.f[creep_id] .= CuArray(f_creep_cpu)

        # Step 2: Compute K*f (boundary displacement contribution to RHS)
        apply_stiffness!(state.a, state.f, weights_gpu, H_gpu, Ht_gpu, dof_id_gpu, mesh.n_elements)
        rhs_full = state.a   # K * f_prescribed
        rhs_gpu = Array(rhs_full[fltni_gpu])  # reduced RHS on CPU for AMGX upload

        # Warm-start from previous displacement
        x0_gpu = Array(state.u_prev[fltni_gpu])

        # Step 2b: AMGX solve: K*u = -rhs (PCG+AMG fully on GPU)
        AMGX.upload!(solver.b_buf, -rhs_gpu)
        AMGX.upload!(solver.x_buf, x0_gpu)
        AMGX.solve!(solver.x_buf, solver.amgx_solver, solver.b_buf)

        # Extract solution back to GPU state vector
        u_sol_cpu = AMGX.download(solver.x_buf)
        if any(isnan.(u_sol_cpu)) || any(isinf.(u_sol_cpu))
            @error "Non-finite values in AMGX solution" pass iteration=state.iteration
            error("AMGX solver produced non-finite values")
        end
        state.u[fltni_gpu] .= CuArray(u_sol_cpu)

        # Enforce prescribed displacements on boundaries (scatter-upload)
        state.u[fault_id] .= CuArray(f_fault_cpu)
        state.u[creep_id] .= CuArray(f_creep_cpu)

        # Step 3: Compute fault traction
        fill!(state.a, zero(T))
        apply_stiffness!(state.a, state.u, weights_gpu, H_gpu, Ht_gpu, dof_id_gpu, mesh.n_elements)
        state.a[creep_id] .= 0

        # Extract fault stress (download fault-size arrays to CPU for NR)
        a_fault_cpu = Array(state.a[fault_id])
        τf_cpu = -a_fault_cpu ./ fault_matrix

        if any(isnan.(τf_cpu)) || any(isinf.(τf_cpu))
            @error "Non-finite fault stress" pass iteration=state.iteration
            error("Non-finite fault stress computed")
        end

        # Steps 4 & 5: Rate-state evolution (CPU loop over ~1-2k fault nodes)
        ψ_cpu = Array(state.ψ)

        for i in eachindex(fault_id)
            ψ_cpu[i] = state_time_evolution(ψ_cpu[i], Vf_new[i], dt,
                                            ics.friction.Lc[i], params.Vo)
            Vf_new[i] = fault_slip_rate(ψ_cpu[i], τf_cpu[i],
                                        ics.τo[i], ics.σo[i],
                                        ics.friction.a[i], ics.friction.b[i],
                                        params.Vo, params.fo)

            if !isfinite(Vf_new[i]) && state.iteration == 1
                @error "Non-finite slip rate at fault node" i pass
                error("Non-finite slip rate computed")
            end
        end

        # Upload updated ψ, τf back to GPU
        copyto!(state.ψ, CuArray(ψ_cpu))
        copyto!(state.τf, CuArray(τf_cpu))

        # Average slip rates
        Vf_new .= 0.5 .* (Vf_old .+ Vf_new)

        # Update fault/creep velocity (scatter-upload)
        v_fault_new_cpu = 0.5 .* (Vf_new .- params.Vpl)
        state.v[fault_id] .= CuArray(v_fault_new_cpu)
        state.v[creep_id] .= zero(T)
    end  # Two-pass iteration

    # Update global velocity from displacement change
    state.v .= (state.u .- state.u_prev) ./ dt
    clamp!(state.v, -physics.vs, physics.vs)

    # Re-enforce boundary velocities (scatter-upload)
    v_fault_final = 0.5 .* (Vf_new .- params.Vpl)
    state.v[fault_id] .= CuArray(v_fault_final)
    state.v[creep_id] .= zero(T)

    # Blowup detection (download max for check)
    v_max = maximum(abs.(Array(state.v)))
    u_max = maximum(abs.(Array(state.u)))
    if v_max > 1e10 || u_max > 1e10
        @error "GPU QS blowup detected" v_max u_max dt iteration=state.iteration
        error("GPU quasistatic solver produced non-physical values")
    end

    fill!(state.a, zero(T))
    state.u[creep_id] .= zero(T)
    state.v[creep_id] .= zero(T)

    copyto!(state.Vf, CuArray(Vf_new))

    return nothing
end
