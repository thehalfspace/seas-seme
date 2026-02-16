"""
    FaultSolvers

Rate-state friction fault boundary solvers for SEAS simulations.

This module provides:
- Slip rate calculation from rate-state friction law
- State variable time evolution (aging law)
- Newton-Raphson search for dynamic fault boundary conditions

Based on Kaneko et al. (2008, 2011) formulation.
"""

"""
    fault_slip_rate(ψ, τ, τ_init, σ_n, a, b, V₀, f₀)

Compute fault slip rate from rate-state friction law (quasi-static formulation).

# Arguments
- `ψ`: Transformed state variable (log θ)
- `τ`: Current shear stress on fault
- `τ_init`: Initial shear stress
- `σ_n`: Normal stress (positive in compression)
- `a`: Direct effect parameter
- `b`: Evolution effect parameter
- `V₀`: Reference slip rate
- `f₀`: Reference friction coefficient

# Returns
- Fault slip rate `V`

# Notes
Uses the stick-traction formulation from Kaneko et al. (2008, 2011).
The `fact` parameter (velocity weakening factor) is set to 1 for quasi-static regime.
"""
function fault_slip_rate(ψ::T, τ::T, τ_init::T, σ_n::T,
                        a::T, b::T, V₀::T, f₀::T) where T<:Real
    # Safety checks for numerical stability
    if a <= 0
        error("Rate-state parameter 'a' must be positive, got a=$a")
    end
    if σ_n <= 0
        error("Normal stress must be positive (compressive), got σ_n=$σ_n")
    end

    fact = one(T)  # For quasi-static: fact = 1

    term1 = fact * (τ + τ_init) / (σ_n * a)
    term2 = -(f₀ + b * ψ) / a

    # Clamp exponential arguments to prevent overflow
    # exp(700) ≈ 1e304, exp(-700) ≈ 3e-305 for Float64
    max_exp_arg = T(700)
    arg_plus = clamp(term2 + term1, -max_exp_arg, max_exp_arg)
    arg_minus = clamp(term2 - term1, -max_exp_arg, max_exp_arg)

    return V₀ * (exp(arg_plus) - exp(arg_minus))
end


"""
    state_time_evolution(ψ_old, V, dt, Lc, V₀)

Time evolution of state variable using aging law (IDState = 2 in Kaneko's code).

# Arguments
- `ψ_old`: State variable at previous time
- `V`: Slip rate
- `dt`: Time step
- `Lc`: Critical slip distance
- `V₀`: Reference slip rate

# Returns
- Updated state variable `ψ_new`

# Notes
This uses a specialized time-integration scheme from Kaneko's Matlab code,
different from standard aging law integration. Handles both small and large
slip regimes with appropriate approximations.
"""
function state_time_evolution(ψ_old::T, V::T, dt::T, Lc::T, V₀::T) where T<:Real
    # Safety checks
    if Lc <= 0
        error("Critical slip distance must be positive, got Lc=$Lc")
    end
    if dt < 0
        error("Time step must be non-negative, got dt=$dt")
    end

    # Clamp slip rate to physical bounds before state evolution
    # |V| > vs (~3500 m/s) is unphysical; use a generous bound
    V_clamped = clamp(V, T(-1e4), T(1e4))

    Vdt_L = abs(V_clamped) * dt / Lc

    if Vdt_L < V₀
        # Small slip regime
        # Clamp exponential argument to prevent overflow
        exp_arg = clamp(ψ_old - Vdt_L, T(-700), T(700))
        log_arg = exp(exp_arg) + V₀ * dt / Lc - 0.5 * V₀ * abs(V_clamped) * dt^2 / Lc^2

        # Ensure log argument is positive
        if log_arg <= 0
            @warn "Non-positive argument to log in state evolution (small slip)" log_arg V dt Lc
            return ψ_old  # Return previous state as fallback
        end

        # Clamp ψ to physical bounds: ψ ∈ [-20, 30]
        # ψ = -20 → V/V₀ ~ e^20 ~ 5e8 (way beyond seismic)
        # ψ = +30 → V/V₀ ~ e^-30 ~ 1e-13 (fully locked)
        return clamp(log(log_arg), T(-20), T(30))
    else
        # Large slip regime
        exp_arg1 = clamp(ψ_old - Vdt_L, T(-700), T(700))
        exp_arg2 = clamp(-Vdt_L, T(-700), T(700))

        log_arg = exp(exp_arg1) + (V₀ / abs(V_clamped)) * (1 - exp(exp_arg2))

        # Ensure log argument is positive
        if log_arg <= 0
            @warn "Non-positive argument to log in state evolution (large slip)" log_arg V dt Lc
            return ψ_old  # Return previous state as fallback
        end

        return clamp(log(log_arg), T(-20), T(30))
    end
end


"""
    nr_search(τ_guess, f₀, V₀, a, b, σ_n, τ_init, ψ, fault_z, fault_vfree;
              tol_factor=0.001, max_iter=1000)

Newton-Raphson search for dynamic fault boundary conditions.

Solves the coupled system:
1. Stick traction: `fault_z * fault_vfree = fault_z * V + τ - τ_init`
2. Rate-state friction law: `V = f(τ, ψ, ...)`

# Arguments
- `τ_guess`: Initial guess for shear stress
- `f₀`: Reference friction coefficient
- `V₀`: Reference slip rate
- `a`, `b`: Rate-state parameters
- `σ_n`: Normal stress
- `τ_init`: Initial shear stress
- `ψ`: Current state variable
- `fault_z`: Fault impedance
- `fault_vfree`: Stick traction (free velocity)
- `tol_factor`: Tolerance factor (default: 0.001)
- `max_iter`: Maximum iterations (default: 1000)

# Returns
- `(V, τ)`: Slip rate and shear stress satisfying both equations

# Notes
Convergence tolerance is `tol_factor * a * σ_n`. If convergence fails,
returns current values with error message.
"""
function nr_search(τ_guess::T, f₀::T, V₀::T, a::T, b::T, σ_n::T,
                  τ_init::T, ψ::T, fault_z::T, fault_vfree::T;
                  tol_factor::T=0.001, max_iter::Int=1000) where T<:Real

    # Velocity weakening factor (large Vw ensures fact ≈ 1)
    Vw = T(1.0e10)
    fact = 1 + (V₀ / Vw) * exp(-ψ)

    # Convergence tolerance
    tol = tol_factor * a * σ_n

    # Maximum safe exponential argument (exp(700) ≈ 1e304 for Float64)
    max_exp_arg = T(700)

    # Initialize
    τ = τ_guess
    δ = typemax(T)  # Large initial value
    iter = 0

    # Newton-Raphson iteration
    while abs(δ) > tol
        # Compute slip rate and its derivative
        fa = fact * τ / (σ_n * a)
        help = -(f₀ + b * ψ) / a

        # CRITICAL FIX: Clamp exponential arguments to prevent overflow → NaN
        arg_plus = clamp(help + fa, -max_exp_arg, max_exp_arg)
        arg_minus = clamp(help - fa, -max_exp_arg, max_exp_arg)

        help1 = exp(arg_plus)
        help2 = exp(arg_minus)

        V = V₀ * (help1 - help2)
        V_prime = fact * (V₀ / (a * σ_n)) * (help1 + help2)

        # Newton-Raphson update
        δ = (fault_z * fault_vfree - fault_z * V + τ_init - τ) /
            (1 + fault_z * V_prime)
        τ += δ

        iter += 1

        # Check for divergence or NaN
        if !isfinite(τ) || !isfinite(δ) || abs(δ) > 1e10 || iter >= max_iter
            @warn "NR search failed to converge" iter δ tol τ V σ_n a b ψ fault_z fault_vfree
            # Return safe fallback: previous slip rate and initial stress
            # instead of diverged garbage values
            return fault_vfree, τ_init
        end
    end

    # Recompute final slip rate with converged stress (with clamping for safety)
    fa = fact * τ / (σ_n * a)
    help = -(f₀ + b * ψ) / a
    arg_plus = clamp(help + fa, -max_exp_arg, max_exp_arg)
    arg_minus = clamp(help - fa, -max_exp_arg, max_exp_arg)
    help1 = exp(arg_plus)
    help2 = exp(arg_minus)
    V = V₀ * (help1 - help2)

    return V, τ
end
