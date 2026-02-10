"""
    FaultSolversPlaneStrain

Fault boundary solvers for plane-strain formulation.

The rate-state friction law is scalar (shear traction vs slip rate), so the
core NR solver from FaultSolvers.jl is reused. This module provides the
projection layer that decomposes 2-component forces/velocities into
fault-tangential and fault-normal components.
"""

"""
    fault_traction_from_force_plane_strain(force, fault_id, fault_mat,
                                           tangent, normal, ndof)

Compute fault shear traction from internal forces (plane-strain).

# Arguments
- `force::Vector{T}`: Internal force vector [2*ndof], component-major
- `fault_id::Vector{Int}`: Fault node spatial DOF indices
- `fault_mat::Vector{T}`: Fault boundary impedance matrix
- `tangent::Matrix{T}`: Unit tangent vectors [2, nfault]
- `normal::Matrix{T}`: Unit normal vectors [2, nfault]
- `ndof::Int`: Number of spatial DOFs

# Returns
- `τ_shear::Vector{T}`: Shear traction (tangential) at each fault node
"""
function fault_traction_from_force_plane_strain(
    force::Vector{T},
    fault_id::Vector{Int},
    fault_mat::Vector{T},
    tangent::Matrix{T},
    normal::Matrix{T},
    ndof::Int
) where {T<:AbstractFloat}
    nfault = length(fault_id)
    τ_shear = zeros(T, nfault)

    for i in 1:nfault
        nid = fault_id[i]
        fx = force[nid]
        fy = force[ndof + nid]

        # Project force onto tangent direction and divide by impedance
        τ_shear[i] = -(fx * tangent[1, i] + fy * tangent[2, i]) / fault_mat[i]
    end

    return τ_shear
end


"""
    compute_stick_traction_plane_strain(v, a, M_global, fault_id,
                                        tangent, ndof, dt)

Compute tangential stick traction (free velocity) for dynamic solver.

FaultVFree_tangential = 2*v_tang + dt*a_tang/M

# Arguments
- `v::Vector{T}`: Velocity field [2*ndof]
- `a::Vector{T}`: Force/acceleration field [2*ndof]
- `M_global::Vector{T}`: Mass matrix [2*ndof]
- `fault_id::Vector{Int}`: Fault node spatial DOF indices
- `tangent::Matrix{T}`: Unit tangent vectors [2, nfault]
- `ndof::Int`: Number of spatial DOFs
- `dt::T`: Timestep

# Returns
- `vfree_tang::Vector{T}`: Tangential free velocity at fault nodes
"""
function compute_stick_traction_plane_strain(
    v::Vector{T},
    a::Vector{T},
    M_global::Vector{T},
    fault_id::Vector{Int},
    tangent::Matrix{T},
    ndof::Int,
    dt::T
) where {T<:AbstractFloat}
    nfault = length(fault_id)
    vfree_tang = zeros(T, nfault)

    for i in 1:nfault
        nid = fault_id[i]
        tx, ty = tangent[1, i], tangent[2, i]

        # Velocity tangential component
        v_tang = v[nid] * tx + v[ndof + nid] * ty

        # Acceleration tangential component (force / mass)
        a_tang = (a[nid] * tx) / M_global[nid] + (a[ndof + nid] * ty) / M_global[ndof + nid]

        vfree_tang[i] = 2 * v_tang + dt * a_tang
    end

    return vfree_tang
end


"""
    apply_fault_traction_plane_strain!(a, fault_id, fault_mat, τf,
                                       tangent, ndof)

Apply fault shear traction to force vector (plane-strain).

The traction is applied in the tangent direction:
  a[fault_x] -= fault_mat * τf * tangent_x
  a[fault_y] -= fault_mat * τf * tangent_y

# Arguments
- `a::Vector{T}`: Force/acceleration vector [2*ndof] (modified in-place)
- `fault_id::Vector{Int}`: Fault node spatial DOF indices
- `fault_mat::Vector{T}`: Fault boundary impedance
- `τf::Vector{T}`: Fault perturbation shear traction
- `tangent::Matrix{T}`: Unit tangent vectors [2, nfault]
- `ndof::Int`: Number of spatial DOFs
"""
function apply_fault_traction_plane_strain!(
    a::Vector{T},
    fault_id::Vector{Int},
    fault_mat::Vector{T},
    τf::Vector{T},
    tangent::Matrix{T},
    ndof::Int
) where {T<:AbstractFloat}
    for i in eachindex(fault_id)
        nid = fault_id[i]
        traction_force = fault_mat[i] * τf[i]
        a[nid]        -= traction_force * tangent[1, i]
        a[ndof + nid] -= traction_force * tangent[2, i]
    end
end


"""
    compute_fault_impedance_plane_strain(M_global, fault_id, fault_mat,
                                         ndof, dt)

Compute scalar fault impedance Z for NR search (plane-strain).

For each fault node: Z = M_tang / (fault_mat * dt)
where M_tang is the effective mass in the tangential direction.
For diagonal mass with identical M_x == M_y, this simplifies to M/fault_mat/dt.

# Returns
- `fault_z::Vector{T}`: Fault impedance at each node
"""
function compute_fault_impedance_plane_strain(
    M_global::Vector{T},
    fault_id::Vector{Int},
    fault_mat::Vector{T},
    ndof::Int,
    dt::T
) where {T<:AbstractFloat}
    nfault = length(fault_id)
    fault_z = zeros(T, nfault)

    for i in 1:nfault
        nid = fault_id[i]
        # M_x == M_y for our formulation, use M_x
        fault_z[i] = M_global[nid] / (fault_mat[i] * dt)
    end

    return fault_z
end


"""
    set_fault_velocity_plane_strain!(v, fault_id, Vf, Vpl, tangent, ndof)

Set velocity at fault nodes from scalar slip rate (plane-strain).

v_x[fault] = 0.5*(Vf - Vpl) * tangent_x
v_y[fault] = 0.5*(Vf - Vpl) * tangent_y
"""
function set_fault_velocity_plane_strain!(
    v::Vector{T},
    fault_id::Vector{Int},
    Vf::Vector{T},
    Vpl::T,
    tangent::Matrix{T},
    ndof::Int
) where {T<:AbstractFloat}
    for i in eachindex(fault_id)
        nid = fault_id[i]
        v_scalar = T(0.5) * (Vf[i] - Vpl)
        v[nid]        = v_scalar * tangent[1, i]
        v[ndof + nid] = v_scalar * tangent[2, i]
    end
end


"""
    get_fault_tangential_velocity(v, fault_id, tangent, ndof, Vpl)

Extract total fault slip rate (tangential) from velocity field.

Returns Vf = 2*v_tangential + Vpl for each fault node.
"""
function get_fault_tangential_velocity(
    v::Vector{T},
    fault_id::Vector{Int},
    tangent::Matrix{T},
    ndof::Int,
    Vpl::T
) where {T<:AbstractFloat}
    nfault = length(fault_id)
    Vf = zeros(T, nfault)

    for i in 1:nfault
        nid = fault_id[i]
        v_tang = v[nid] * tangent[1, i] + v[ndof + nid] * tangent[2, i]
        Vf[i] = 2 * v_tang + Vpl
    end

    return Vf
end
