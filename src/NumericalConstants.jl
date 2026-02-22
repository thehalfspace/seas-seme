"""
    NumericalConstants

Named constants for numerical stability in rate-state friction solvers.

All constants are dimensionless or carry well-defined physical meaning.
Using named constants instead of magic numbers makes the bounds auditable
against the literature and adjustable for non-Float64 precision.

# Exponential bounds
- `EXP_ARG_MAX`: Maximum safe argument to `exp()` for Float64.
  `exp(709)` is finite; `exp(710)` overflows to Inf. Using 700 gives
  a 9-decade safety margin. Source: `log(floatmax(Float64)) ≈ 709.78`.

# State variable bounds (ψ = log(θ V₀ / Lc))
- `PSI_MIN`: ψ = -20 corresponds to V/V₀ ~ e²⁰ ~ 5×10⁸ m/s, far beyond
  any physical slip rate. Acts as a floor for extreme velocity events.
- `PSI_MAX`: ψ = +30 corresponds to V/V₀ ~ e⁻³⁰ ~ 1×10⁻¹³, fully locked.
  Ruina (1983); Dieterich (1979).

# Slip rate bounds
- `V_PHYS_MAX`: Generous upper bound on physical slip rate (m/s).
  Seismic rupture velocities reach ~5 m/s; 1e4 m/s ≈ 3×Vp and is
  far beyond anything physical. Used to clamp V before state evolution.

# Weakening velocity
- `VW_WEAKENING`: Reference velocity for the flash-heating weakening factor
  `fact = 1 + (V₀/Vw) * exp(-ψ)`. Setting Vw = 1e10 makes fact ≈ 1
  for all physical slip rates, effectively disabling weakening while
  preserving the formula's structure for future extension.
"""

# Maximum argument to exp() that does not overflow Float64
const EXP_ARG_MAX = 700.0

# State variable ψ clamp bounds (Dieterich 1979, Ruina 1983)
const PSI_MIN = -20.0
const PSI_MAX =  30.0

# Physical slip rate upper bound for clamping before state evolution (m/s)
const V_PHYS_MAX = 1e4

# Flash-heating weakening velocity — set large so fact ≈ 1 (weakening disabled)
const VW_WEAKENING = 1e10
