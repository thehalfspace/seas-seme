"""
    DynamicSolverPlaneStrain

Explicit dynamic earthquake solver for plane-strain formulation.
Leap-frog time integration on 2*ndof displacement/velocity fields.

Based on Kaneko et al. (2008) Appendix B, extended to 2-component in-plane motion.

Phase 1: Uses S-wave absorbing BC (simplified).
Phase 3 will upgrade to P+SV directional absorbing BC.
"""

"""
    DynamicSolverPlaneStrain{T}

Explicit dynamic solver for plane-strain.

# Fields
- `dt_min::T`: Minimum timestep from CFL condition [s]
- `verbose::Bool`: Print diagnostic information

# Notes
For plane-strain, CFL should be based on P-wave velocity (faster than S-wave):
  dt_min = CFL * dx / Vp
This is handled in build_simulation_plane_strain.
"""
struct DynamicSolverPlaneStrain{T<:AbstractFloat}
    dt_min::T
    verbose::Bool
end

function DynamicSolverPlaneStrain(dt_min::T; verbose::Bool=false) where T<:AbstractFloat
    return DynamicSolverPlaneStrain{T}(dt_min, verbose)
end


"""
    dynamic_step!(state::SimulationStatePlaneStrain, solver, mesh, physics, ics,
                  params, M_global, weights, H, Ht, dof_id)

Perform one dynamic time step using leap-frog integration (plane-strain).

# Algorithm (Kaneko et al. 2008, Appendix B, extended to 2-component)
1. Update displacement: u[n+1] = u[n] + dt*v[n] + 0.5*dt²*a[n]
2. Partial velocity: v[n+1/2] = v[n] + 0.5*dt*a[n]
3. Internal forces: f = -K*u[n+1] (via tensor-product matvec)
4. Absorbing BC (S-wave damping, both components)
5. Compute stick traction (tangential free velocity)
6. Update state variable
7. Newton-Raphson search (two iterations)
8. Apply fault traction
9. Compute acceleration: a = f / M
10. Complete velocity: v[n+1] = v[n+1/2] + 0.5*dt*a[n+1]
"""
function dynamic_step!(state::SimulationStatePlaneStrain{T},
                       solver::DynamicSolverPlaneStrain{T},
                       mesh, physics, ics, params,
                       M_global::AbstractVector{T},
                       weights::MetricWeightsPlaneStrain{T},
                       H::AbstractMatrix{T}, Ht::AbstractMatrix{T},
                       dof_id::AbstractArray{<:Integer,3}) where T<:AbstractFloat

    dt = solver.dt_min
    ndof = state.ndof

    fault_id = mesh.boundaries.fault.node_ids
    creep_id = mesh.boundaries.creep.node_ids
    absorbing_id = mesh.boundaries.absorbing.node_ids
    fault_matrix = mesh.boundaries.fault.matrix
    mask = mesh.active_fault_mask
    tangent = state.fault_tangent
    normal = state.fault_normal

    absorb_matrix = mesh.boundaries.absorbing.matrix
    absorb_tangent = mesh.boundaries.absorbing.tangent
    absorb_normal = mesh.boundaries.absorbing.normal

    # Fault impedance
    fault_z = compute_fault_impedance_plane_strain(M_global, fault_id, fault_matrix, ndof, dt)

    # Store previous state
    copyto!(state.u_prev, state.u)
    copyto!(state.v_prev, state.v)

    # Step 1: Update displacement (leap-frog)
    state.u .= state.u .+ dt .* state.v .+ (T(0.5) * dt^2) .* state.a

    # Step 2: Partial velocity update
    state.v .= state.v .+ (T(0.5) * dt) .* state.a
    fill!(state.a, zero(T))

    # Step 3: Internal forces -K*u via tensor-product matvec
    apply_stiffness_plane_strain!(state.a, state.u, weights, H, Ht, dof_id, mesh.n_elements, ndof)
    state.a .= -state.a

    # Zero force on creep boundary
    state.a[creep_id] .= 0
    state.a[ndof .+ creep_id] .= 0

    # Step 4: Absorbing BC (S-wave damping for both components)
    # Phase 1: simplified — treat as scalar impedance per component
    # For GPU: download absorbing nodes, modify, upload (absorbing_id is small)
    if !(state.a isa Vector)
        a_abs_x = Array(state.a[absorbing_id])
        a_abs_y = Array(state.a[ndof .+ absorbing_id])
        v_abs_x = Array(state.v[absorbing_id])
        v_abs_y = Array(state.v[ndof .+ absorbing_id])
        for i in eachindex(absorbing_id)
            a_abs_x[i] -= absorb_matrix[i] * v_abs_x[i]
            a_abs_y[i] -= absorb_matrix[i] * v_abs_y[i]
        end
        state.a[absorbing_id]        .= CuArray(a_abs_x)
        state.a[ndof .+ absorbing_id] .= CuArray(a_abs_y)
    else
        @inbounds for i in eachindex(absorbing_id)
            nid = absorbing_id[i]
            state.a[nid]        -= absorb_matrix[i] * state.v[nid]
            state.a[ndof + nid] -= absorb_matrix[i] * state.v[ndof + nid]
        end
    end

    # Step 5: Compute stick traction (tangential free velocity)
    # compute_stick_traction_plane_strain downloads from GPU internally
    fault_vfree_cpu = compute_stick_traction_plane_strain(
        state.v, state.a, M_global, fault_id, tangent, ndof, dt
    )

    # Previous fault slip rate (tangential) — downloads GPU arrays internally
    Vf_prev = get_fault_tangential_velocity(state.v_prev, fault_id, tangent, ndof, params.Vpl)

    # Download fault-size arrays to CPU for NR loop
    ψ_cpu   = Array(state.ψ)
    Vf_cpu  = Array(state.Vf)
    τf_cpu  = Array(state.τf)

    # Workspace for two-iteration NR
    Vf_first_iter = similar(Vf_cpu)

    # Steps 5-7: Solve nonlinear fault boundary equations (CPU loop)
    for i in eachindex(fault_id)
        if !mask[i]
            # Excluded endpoint: lock to plate rate, no state evolution
            Vf_cpu[i] = params.Vpl
            τf_cpu[i] = ics.τo[i]  # total stress = initial (no perturbation)
            Vf_first_iter[i] = params.Vpl
            continue
        end

        # Normal stress: use constant σ_n⁰ for now.
        # TODO: For dipping faults, compute Δσ_n from elastic forces and use
        # σn_eff = ics.σo[i] + state.σn_perturbation[i]
        σn_eff = ics.σo[i]

        # Save ψ before update so we can restore on NR failure
        ψ_saved = ψ_cpu[i]

        # Update state variable
        ψ_cpu[i] = state_time_evolution(ψ_cpu[i], Vf_prev[i], dt,
                                        ics.friction.Lc[i], params.Vo)

        Vf_first_iter[i] = Vf_prev[i]

        # Convert perturbation stress to total for NR
        τ_total_guess = τf_cpu[i] + ics.τo[i]

        # Newton-Raphson search (1st iteration)
        Vf_cpu[i], τf_cpu[i] = nr_search(
            τ_total_guess, params.fo, params.Vo,
            ics.friction.a[i], ics.friction.b[i], σn_eff,
            ics.τo[i], ψ_cpu[i], fault_z[i], fault_vfree_cpu[i]
        )

        if isnan(Vf_cpu[i]) || isnan(τf_cpu[i])
            # NR failed completely — restore ψ and use previous slip rate
            ψ_cpu[i] = ψ_saved
            Vf_cpu[i] = Vf_prev[i]
            τf_cpu[i] = τ_total_guess - ics.τo[i]
            @warn "NR search produced NaN (plane-strain), restoring previous state" location=i iter=state.iteration
            continue
        end

        # 2nd iteration with average slip rate
        ψ_after_first = ψ_cpu[i]
        avg_slip_rate = T(0.5) * (Vf_cpu[i] + Vf_first_iter[i])
        ψ_cpu[i] = state_time_evolution(ψ_cpu[i], avg_slip_rate, dt,
                                        ics.friction.Lc[i], params.Vo)

        Vf_2nd, τf_2nd = nr_search(
            τf_cpu[i], params.fo, params.Vo,
            ics.friction.a[i], ics.friction.b[i], σn_eff,
            ics.τo[i], ψ_cpu[i], fault_z[i], fault_vfree_cpu[i]
        )

        if isnan(Vf_2nd) || isnan(τf_2nd)
            # 2nd NR failed — keep 1st iteration results and restore ψ
            ψ_cpu[i] = ψ_after_first
            @warn "NR 2nd iteration produced NaN, keeping 1st iteration result" location=i
        else
            Vf_cpu[i] = Vf_2nd
            τf_cpu[i] = τf_2nd
        end
    end

    # Upload NR results back to state (CPU or GPU)
    copyto!(state.ψ, _as_storage(state.ψ, ψ_cpu))
    copyto!(state.Vf, _as_storage(state.Vf, Vf_cpu))
    copyto!(state.fault_vfree, _as_storage(state.fault_vfree, fault_vfree_cpu))

    # Convert to perturbation stress and upload
    τf_pert_cpu = τf_cpu .- ics.τo
    copyto!(state.τf, _as_storage(state.τf, τf_pert_cpu))

    # Step 8: Apply fault traction to force vector (tangential direction)
    # apply_fault_traction_plane_strain! handles GPU arrays internally
    apply_fault_traction_plane_strain!(state.a, fault_id, fault_matrix, state.τf,
                                       tangent, ndof)

    # Step 9: Compute acceleration
    state.a .= state.a ./ M_global

    # Step 10: Complete velocity update
    state.v .= state.v .+ (T(0.5) * dt) .* state.a

    # Enforce zero velocity on creep boundary
    state.v[creep_id] .= 0
    state.v[ndof .+ creep_id] .= 0
    state.a[creep_id] .= 0
    state.a[ndof .+ creep_id] .= 0

    return nothing
end
