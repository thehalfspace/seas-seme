"""
    SimulationStatePlaneStrain

Container for all simulation state variables for plane-strain formulation.

Global fields use 2*ndof vectors in component-major order: [u_x(1..ndof), u_y(1..ndof)].
Fault fields remain scalar (tangential slip rate, shear traction) since the
rate-state friction law operates on scalar traction and slip rate.
"""

"""
    SimulationStatePlaneStrain{T<:AbstractFloat}

Mutable container for plane-strain simulation state.

# Global Fields (length 2*ndof: [component_x; component_y])
- `u::Vector{T}`: Displacement field [u_x; u_y]
- `v::Vector{T}`: Velocity field [v_x; v_y]
- `a::Vector{T}`: Acceleration/force field [a_x; a_y]

# Fault-Specific Fields (length nfault, scalar tangential quantities)
- `τf::Vector{T}`: Fault shear traction (tangential component)
- `ψ::Vector{T}`: Transformed state variable (log θ)
- `Vf::Vector{T}`: Fault slip rate (tangential, scalar)
- `σn_perturbation::Vector{T}`: Normal stress perturbation from elasticity (Δσ_n)

# Workspace Arrays
- `u_prev::Vector{T}`: Previous displacement (2*ndof)
- `v_prev::Vector{T}`: Previous velocity (2*ndof)
- `f::Vector{T}`: RHS vector for linear solver (2*ndof)
- `fault_vfree::Vector{T}`: Stick traction, tangential (nfault)

# Fault Geometry (precomputed, constant)
- `fault_tangent::Matrix{T}`: Unit tangent vectors [2, nfault]
- `fault_normal::Matrix{T}`: Unit outward normal vectors [2, nfault]

# Time Tracking
- `time::T`: Current simulation time [s]
- `timestep::T`: Current timestep [s]
- `iteration::Int`: Iteration counter
- `solver_mode::Symbol`: Current solver (:quasistatic or :dynamic)

# Dimensions
- `ndof::Int`: Number of spatial DOFs
- `nfault::Int`: Number of fault nodes
"""
mutable struct SimulationStatePlaneStrain{T<:AbstractFloat}
    # Global fields (length 2*ndof)
    u::Vector{T}
    v::Vector{T}
    a::Vector{T}

    # Fault-specific fields (length nfault)
    τf::Vector{T}
    ψ::Vector{T}
    Vf::Vector{T}
    σn_perturbation::Vector{T}  # Normal stress perturbation from elasticity

    # Workspace arrays
    u_prev::Vector{T}
    v_prev::Vector{T}
    f::Vector{T}
    fault_vfree::Vector{T}

    # Fault geometry
    fault_tangent::Matrix{T}    # [2, nfault]
    fault_normal::Matrix{T}     # [2, nfault]

    # Time tracking
    time::T
    timestep::T
    iteration::Int
    solver_mode::Symbol

    # Dimensions
    ndof::Int
    nfault::Int
end


"""
    SimulationStatePlaneStrain(mesh, ics, params; v_init=5.0e-4)

Initialize plane-strain simulation state.

# Arguments
- `mesh::UnstructuredSEMesh`: Spectral element mesh
- `ics`: Initial conditions (stress, friction parameters)
- `params`: Simulation parameters (Vpl, V₀, f₀)
- `v_init::Real`: Initial velocity magnitude [m/s] (default: 0.5 mm/s)

# Notes
Initial velocity is projected onto the fault tangent direction.
For a vertical fault (dip=90°), tangent = (0, -1), so the initial velocity
is entirely in the -y direction, consistent with downward-sliding motion.
"""
function SimulationStatePlaneStrain(mesh::UnstructuredSEMesh{T}, ics, params;
                                    v_init::Real=T(5.0e-4)) where T<:AbstractFloat
    ndof = mesh.ndof
    nfault = length(mesh.boundaries.fault.node_ids)
    fault_id = mesh.boundaries.fault.node_ids
    creep_id = mesh.boundaries.creep.node_ids

    # Get fault tangent and normal from mesh boundary data
    fault_tangent = copy(mesh.boundaries.fault.tangent)   # [2, nfault]
    fault_normal = copy(mesh.boundaries.fault.normal)     # [2, nfault]

    # Allocate global fields (2*ndof)
    u = zeros(T, 2 * ndof)
    v = zeros(T, 2 * ndof)
    a = zeros(T, 2 * ndof)

    # Allocate fault fields
    τf = zeros(T, nfault)
    ψ = zeros(T, nfault)
    Vf = zeros(T, nfault)
    σn_perturbation = zeros(T, nfault)

    # Allocate workspace arrays
    u_prev = zeros(T, 2 * ndof)
    v_prev = zeros(T, 2 * ndof)
    f = zeros(T, 2 * ndof)
    fault_vfree = zeros(T, nfault)

    # Initialize velocity field
    # For single-sided fault: v = 0.5*(Vf - Vpl) projected onto tangent
    # Start with uniform initial slip rate v_init
    v_fault_scalar = T(v_init) - T(0.5) * params.Vpl

    # Set velocity at fault nodes (projected onto tangent direction)
    for i in 1:nfault
        nid = fault_id[i]
        v[nid]        = v_fault_scalar * fault_tangent[1, i]   # v_x
        v[ndof + nid] = v_fault_scalar * fault_tangent[2, i]   # v_y
    end

    # Creep boundary: zero velocity
    v[creep_id] .= 0
    v[ndof .+ creep_id] .= 0

    # Initialize fault slip rate (scalar tangential)
    # Vf = 2*v_tangential + Vpl
    for i in 1:nfault
        nid = fault_id[i]
        v_tang = v[nid] * fault_tangent[1, i] + v[ndof + nid] * fault_tangent[2, i]
        Vf[i] = 2 * v_tang + params.Vpl
    end

    # Initialize transformed state variable ψ from BP3-FD eq. 23.
    # ψ = (a/b) * ln( (2*V₀/V_init) * sinh(τ⁰/(a*σ_n)) ) - f₀/b
    # This ensures ψ is consistent with τ⁰ and the regularized friction law.
    # V_init = Vf (initial slip rate) for each node.
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

    # Enforce excluded fault nodes to plate rate (Vf = Vpl → v = 0)
    mask = mesh.active_fault_mask
    for i in 1:nfault
        if !mask[i]
            Vf[i] = params.Vpl
            nid = fault_id[i]
            v[nid]        = zero(T)
            v[ndof + nid] = zero(T)
        end
    end

    # Time tracking
    time = zero(T)
    timestep = zero(T)
    iteration = 0
    solver_mode = :quasistatic

    return SimulationStatePlaneStrain{T}(
        u, v, a,
        τf, ψ, Vf, σn_perturbation,
        u_prev, v_prev, f, fault_vfree,
        fault_tangent, fault_normal,
        time, timestep, iteration, solver_mode,
        ndof, nfault
    )
end


"""
    maximum_fault_slip_rate(state::SimulationStatePlaneStrain)

Get maximum fault slip rate from current state.
"""
function maximum_fault_slip_rate(state::SimulationStatePlaneStrain)
    return maximum(abs.(state.Vf))
end
