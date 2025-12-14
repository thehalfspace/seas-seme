"""
Rate-state friction law types and parameters.

Implements the aging law formulation of rate-and-state friction
(Dieterich 1979, Ruina 1983).
"""

"""
    RateStateFriction{T}

Rate-and-state friction parameters with aging law.

# Fields
- `a::Vector{T}`: Direct effect parameter (dimensionless)
- `b::Vector{T}`: Evolution effect parameter (dimensionless)
- `Lc::Vector{T}`: Critical slip distance (m)
- `Vo::T`: Reference slip rate (m/s)
- `fo::T`: Reference friction coefficient (dimensionless)

# Friction Law
The friction coefficient μ depends on slip rate V and state variable θ:
```
μ = a * asinh(V/(2*Vo) * exp(θ*Vo/Lc))
```

Or equivalently with transformed state variable ψ = log(θ*Vo/Lc):
```
μ = a * asinh((V/(2*Vo)) * exp(ψ))
```

# State Evolution (Aging Law)
```
dθ/dt = 1 - V*θ/Lc
```

# Friction Regimes
- Velocity weakening (VW): b > a (seismogenic, unstable)
- Velocity strengthening (VS): a > b (aseismic, stable creep)

# Notes
- Parameters (a, b, Lc) are typically depth-dependent
- Reference values: fo ≈ 0.6, Vo ≈ 10⁻⁶ m/s
"""
struct RateStateFriction{T<:AbstractFloat}
    a::Vector{T}
    b::Vector{T}
    Lc::Vector{T}
    Vo::T
    fo::T
end

"""
    MaterialProperties{T}

Material properties for elastic wave propagation.

# Fields
- `ρ::T`: Density (kg/m³)
- `vs::T`: Shear wave velocity (m/s)
- `μ::T`: Shear modulus (Pa), computed as ρ*vs²

# Notes
Typical crustal values:
- ρ ≈ 2670 kg/m³
- vs ≈ 3464 m/s
- μ ≈ 32 GPa
"""
struct MaterialProperties{T<:AbstractFloat}
    ρ::T
    vs::T
    μ::T
end

"""
    MaterialProperties(ρ::T, vs::T) where T -> MaterialProperties{T}

Construct material properties, computing shear modulus μ = ρ*vs².
"""
function MaterialProperties(ρ::T, vs::T) where {T<:AbstractFloat}
    μ = ρ * vs^2
    return MaterialProperties(ρ, vs, μ)
end

"""
    InitialConditions{T}

Initial conditions for fault simulation.

# Fields
- `σn::Vector{T}`: Initial normal stress on fault (Pa)
- `τs::Vector{T}`: Initial shear stress on fault (Pa)
- `friction::RateStateFriction{T}`: Rate-state friction parameters

# Notes
- Stresses are positive in compression
- Initial stresses set the pre-stress state before simulation
- Combined with friction parameters, determines nucleation behavior
"""
struct InitialConditions{T<:AbstractFloat}
    σn::Vector{T}
    τs::Vector{T}
    friction::RateStateFriction{T}
end
