"""
    Timestepping

Adaptive timestep computation for earthquake cycle simulations.

Implements the algorithm from Lapusta et al. (2000) and Kaneko et al. (2011) Appendix C.
Automatically adjusts timestep based on nucleation length and fault slip velocity.
"""

"""
    AdaptiveTimestepper{T}

Container for adaptive timestepping parameters.

# Fields
- `dt_min::T`: Minimum timestep (CFL condition for dynamic solver)
- `dt_max::T`: Maximum timestep (user-specified limit)
- `hcell::T`: Characteristic cell size (minimum element dimension)
- `μ_max::T`: Maximum shear modulus
- `Vthreshold_qs_to_dyn::T`: Velocity threshold for quasi-static → dynamic
- `Vthreshold_dyn_to_qs::T`: Velocity threshold for dynamic → quasi-static
"""
struct AdaptiveTimestepper{T<:AbstractFloat}
    dt_min::T
    dt_max::T
    hcell::T
    μ_max::T
    Vthreshold_qs_to_dyn::T
    Vthreshold_dyn_to_qs::T
end


"""
    compute_timestep(fault_V, fault_coords, friction, timestepper, solver_mode,
                    dt_current; ξthf=1.0, ξmax=0.5, dt_incf=1.2)

Compute adaptive timestep for earthquake cycle simulation.

# Arguments
- `fault_V`: Fault slip rate vector [m/s]
- `fault_coords`: Fault node coordinates [2 × nfault]
- `friction`: RateStateFriction object with parameters (a, b, Lc, σ_n)
- `timestepper`: AdaptiveTimestepper parameters
- `solver_mode`: Current solver (`:quasistatic` or `:dynamic`)
- `dt_current`: Current timestep [s]
- `ξthf`: Threshold factor for nucleation length (default: 1.0)
- `ξmax`: Maximum fraction of critical slip (default: 0.5)
- `dt_incf`: Maximum timestep increase factor (default: 1.2)

# Returns
- `dt_new`: New timestep [s]

# Algorithm (Quasi-static mode)
For each fault node, compute:
1. Nucleation length parameter `ξth` based on (a-b)/a and material properties
2. Critical slip distance `ξLf` as fraction of Lc
3. Restrict timestep so slip doesn't exceed critical distance: `dt < ξLf / |V|`
4. Also enforce: `dt ≥ dt_min` and `dt ≤ dt_incf * dt_current` (gradual changes)

# Algorithm (Dynamic mode)
Use minimum timestep from CFL condition: `dt = dt_min`

# Notes
The nucleation length calculation follows Lapusta et al. (2000):
```
r₀ = (μ * Lc) / (a * σ_n) * (1 - (b-a)/a)
```
where the threshold `ξth` ensures stable resolution of the nucleation process.
"""
function compute_timestep(fault_V::AbstractVector{T},
                         fault_coords::AbstractMatrix{T},
                         friction,  # RateStateFriction type
                         timestepper::AdaptiveTimestepper{T},
                         solver_mode::Symbol,
                         dt_current::T;
                         ξthf::T=one(T),
                         ξmax::T=T(0.5),
                         dt_incf::T=T(1.2)) where T<:AbstractFloat

    # Download to CPU if needed (fault array is small: ~1-2k nodes)
    fault_V_cpu = fault_V isa Vector ? fault_V : Array(fault_V)

    if solver_mode == :quasistatic
        dt_new = timestepper.dt_max  # Start with maximum allowed

        # Check each fault node
        for i in eachindex(fault_V_cpu)
            # Compute nucleation length threshold parameter
            # expr1 = -(a - b) / a
            expr1 = -(friction.a[i] - friction.b[i]) / friction.a[i]

            # expr2 = 0.25π * (μ / h) * (Lc / (a * σ_n))
            expr2 = T(0.25) * T(π) * (timestepper.μ_max / timestepper.hcell) *
                   (friction.Lc[i] / (friction.a[i] * friction.σ_n[i]))

            # Nucleation parameter
            r₀ = expr2 - expr1

            # Compute threshold ξth
            discriminant = T(0.25) * r₀^2 - expr2
            if discriminant ≥ 0
                ξth = 1 / r₀
            else
                ξth = 1 - expr1 / expr2
            end

            # Critical slip distance (fraction of Lc)
            if ξthf * ξth > ξmax
                ξLf = ξmax * friction.Lc[i]
            else
                ξLf = ξthf * ξth * friction.Lc[i]
            end

            # Restrict timestep based on slip rate
            if abs(fault_V_cpu[i]) * timestepper.dt_max > ξLf
                dt_cell = ξLf / abs(fault_V_cpu[i])
                if dt_cell < dt_new
                    dt_new = dt_cell
                end
            end
        end

        # Enforce minimum QS timestep (10x CFL) to prevent velocity spike
        # from v = du/dt amplification when dt is tiny at dynamic→QS transition
        dt_min_qs = T(10) * timestepper.dt_min
        if dt_new < dt_min_qs
            dt_new = dt_min_qs
        end

        # Gradual timestep increase (no more than dt_incf * current)
        # Skip this check if dt_current is zero (first iteration)
        if dt_current > zero(T) && dt_new > dt_incf * dt_current
            dt_new = dt_incf * dt_current
        end

        return dt_new

    elseif solver_mode == :dynamic
        # Dynamic mode: use minimum CFL timestep
        return timestepper.dt_min

    else
        error("Unknown solver mode: $solver_mode (expected :quasistatic or :dynamic)")
    end
end


"""
    determine_solver_mode(Vf_max, current_mode, timestepper)

Determine solver mode based on maximum fault slip velocity.

# Arguments
- `Vf_max`: Maximum fault slip rate [m/s]
- `current_mode`: Current solver mode (`:quasistatic` or `:dynamic`)
- `timestepper`: AdaptiveTimestepper with velocity thresholds

# Returns
- New solver mode (`:quasistatic` or `:dynamic`)

# Logic
- QS → Dynamic: if `Vf_max ≥ Vthreshold_qs_to_dyn` (typically 5 mm/s)
- Dynamic → QS: if `Vf_max < Vthreshold_dyn_to_qs` (typically 2 mm/s)
- Otherwise: maintain current mode

Hysteresis prevents rapid mode switching.
"""
function determine_solver_mode(Vf_max::T, current_mode::Symbol,
                              timestepper::AdaptiveTimestepper{T}) where T<:AbstractFloat
    # Mode switching with hysteresis
    if current_mode == :quasistatic && Vf_max < timestepper.Vthreshold_qs_to_dyn ||
       current_mode == :dynamic && Vf_max < timestepper.Vthreshold_dyn_to_qs
        return :quasistatic
    else
        return :dynamic
    end
end
