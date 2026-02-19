"""
    SimulationState

Container for all simulation state variables and workspace arrays.

Designed for zero-allocation time stepping - all arrays are pre-allocated
and reused throughout the simulation.
"""

"""
    SimulationState{T<:AbstractFloat}

Mutable container for simulation state.

# Global Fields (all DOFs)
- `u::Vector{T}`: Displacement field
- `v::Vector{T}`: Velocity field
- `a::Vector{T}`: Acceleration/force field

# Fault-Specific Fields (fault DOFs only)
- `τf::Vector{T}`: Fault shear stress
- `ψ::Vector{T}`: Transformed state variable (log θ)
- `Vf::Vector{T}`: Fault slip rate
- `cum_slip::Vector{T}`: Cumulative fault slip (integrated Vf*dt)

# Workspace Arrays (pre-allocated for zero allocations)
- `u_prev::Vector{T}`: Previous displacement
- `v_prev::Vector{T}`: Previous velocity
- `f::Vector{T}`: RHS vector for linear solver
- `fault_vfree::Vector{T}`: Stick traction (dynamic solver)

# Time Tracking
- `time::T`: Current simulation time [s]
- `timestep::T`: Current timestep [s]
- `iteration::Int`: Iteration counter
- `solver_mode::Symbol`: Current solver (`:quasistatic` or `:dynamic`)

# Design Notes
All arrays are mutable and modified in-place during time stepping to avoid allocations.
The state can be serialized for checkpointing via JLD2.
"""
mutable struct SimulationState{T<:AbstractFloat}
    # Global fields
    u::Vector{T}
    v::Vector{T}
    a::Vector{T}

    # Fault-specific fields
    τf::Vector{T}
    ψ::Vector{T}
    Vf::Vector{T}
    cum_slip::Vector{T}

    # Workspace arrays
    u_prev::Vector{T}
    v_prev::Vector{T}
    f::Vector{T}
    fault_vfree::Vector{T}

    # Time tracking
    time::T
    timestep::T
    iteration::Int
    solver_mode::Symbol
end


"""
    SimulationState(mesh::UnstructuredSEMesh, ics, params; v_init=5.0e-4, T=Float64)

Initialize simulation state with specified initial conditions.

# Arguments
- `mesh::UnstructuredSEMesh`: Spectral element mesh
- `ics`: Initial conditions (stress, friction parameters)
- `params`: Simulation parameters (Vpl, V₀, f₀)
- `v_init::Real`: Initial velocity [m/s] (default: 0.5 mm/s)
- `T::Type`: Floating point type (default: Float64)

# Returns
- `SimulationState{T}`: Initialized state ready for time stepping

# Initialization
- Displacement: zero everywhere
- Velocity: `v_init - 0.5*Vpl` (interseismic loading)
- Fault state variable: Computed from initial stress and velocity
- Solver mode: `:quasistatic`
- Time: 0.0
- Iteration: 0

# Example
```julia
state = SimulationState(mesh, ics, params, v_init=5.0e-4)
```
"""
function SimulationState(mesh::UnstructuredSEMesh{T}, ics, params;
                        v_init::Real=T(5.0e-4)) where T<:AbstractFloat

    ndof = mesh.ndof
    nfault = length(mesh.boundaries.fault.node_ids)
    fault_id = mesh.boundaries.fault.node_ids
    creep_id = mesh.boundaries.creep.node_ids

    # Allocate global fields
    u = zeros(T, ndof)
    v = zeros(T, ndof)
    a = zeros(T, ndof)

    # Allocate fault fields
    τf = zeros(T, nfault)
    ψ = zeros(T, nfault)
    Vf = zeros(T, nfault)
    cum_slip = zeros(T, nfault)

    # Allocate workspace arrays
    u_prev = zeros(T, ndof)
    v_prev = zeros(T, ndof)
    f = zeros(T, ndof)
    fault_vfree = zeros(T, nfault)

    # Initialize velocity field (interseismic loading rate)
    # For single-sided fault: v = 0.5*(Vf - Vpl)
    # Start with uniform initial slip rate v_init
    v .= T(v_init) .- T(0.5) * params.Vpl
    v[creep_id] .= 0

    # Initialize fault slip rate
    Vf .= 2 .* v[fault_id] .+ params.Vpl

    # Initialize transformed state variable ψ from BP3-FD eq. 23.
    # ψ = (a/b) * ln( (2*V₀/V_init) * sinh(τ⁰/(a*σ_n)) ) - f₀/b
    # This uses the regularized friction law and ensures ψ is consistent
    # with τ⁰ at the initial slip rate.
    for i in eachindex(ψ)
        V_init_i = abs(Vf[i])
        V_init_i = max(V_init_i, T(1e-20))  # Guard against zero
        ai = ics.friction.a[i]
        bi = ics.friction.b[i]
        σn_i = ics.σo[i]
        τ0_i = ics.τo[i]
        ψ[i] = (ai / bi) * log((2 * params.Vo / V_init_i) *
               sinh(τ0_i / (ai * σn_i))) - params.fo / bi
    end

    # Time tracking
    time = zero(T)
    timestep = zero(T)
    iteration = 0
    solver_mode = :quasistatic

    return SimulationState{T}(
        u, v, a,
        τf, ψ, Vf, cum_slip,
        u_prev, v_prev, f, fault_vfree,
        time, timestep, iteration, solver_mode
    )
end


"""
    current_fault_slip_velocity(state::SimulationState, fault_id, Vpl)

Compute current total fault slip velocity from state.

# Arguments
- `state::SimulationState`: Current simulation state
- `fault_id::Vector{Int}`: Fault node indices
- `Vpl::Real`: Plate loading velocity

# Returns
- Total fault slip velocity [m/s] (not half-rate)

# Notes
For single-sided fault: `V_total = 2*v + Vpl`
"""
function current_fault_slip_velocity(state::SimulationState, fault_id::Vector{Int}, Vpl::Real)
    return 2 .* state.v[fault_id] .+ Vpl
end


"""
    maximum_fault_slip_rate(state::SimulationState)

Get maximum fault slip rate from current state.

# Returns
- Maximum slip rate [m/s]
"""
function maximum_fault_slip_rate(state::SimulationState)
    return maximum(abs.(state.Vf))
end
