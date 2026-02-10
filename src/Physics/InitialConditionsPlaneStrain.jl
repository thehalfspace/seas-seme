"""
Initial condition setup for plane-strain fault simulations.

Same depth-dependent stress and friction parameters as antiplane,
but depth is computed using the dip angle:
  vertical_depth = along_dip_distance * sin(δ)

For δ=90° (vertical fault), this is identical to the antiplane case.
"""

"""
    build_initial_conditions_plane_strain(config, mesh, dip_angle) -> InitialConditions

Construct initial conditions for plane-strain fault simulation.

# Arguments
- `config::PhysicsConfig`: Physics configuration
- `mesh::UnstructuredSEMesh`: Complete mesh with boundaries
- `dip_angle::Float64`: Fault dip angle in degrees

# Returns
- `InitialConditions{Float64}`: Initial stresses and friction parameters

# Notes
Depth is always measured vertically (not along dip), consistent with
lithostatic stress conventions. For a dipping fault, the y-coordinate
of fault nodes already gives the vertical position, so the depth calculation
is the same as antiplane:
  depth = |fault_y - max(fault_y)|

The dip_angle parameter is carried for future use (e.g., resolving far-field
stress onto the fault plane) but does not change the depth zones.
"""
function build_initial_conditions_plane_strain(
    config::PhysicsConfig,
    mesh::UnstructuredSEMesh{T},
    dip_angle::Float64
) where {T<:AbstractFloat}
    # For Phase 1 (vertical fault), the IC generation is identical to antiplane.
    # The dip_angle is validated but doesn't change the depth-zone logic since
    # depth is always vertical.
    #
    # In Phase 2, for dipping faults, the mesh will have fault nodes at the
    # correct (x,y) positions, and depth = |y - y_surface| is still correct.

    @info "Building plane-strain initial conditions (dip_angle = $(dip_angle)°)"

    # Delegate to the existing depth-dependent IC builder
    # This works because:
    # 1. The depth zones use vertical depth from surface
    # 2. The fault node y-coordinates correctly give vertical position
    # 3. The rate-state parameters are independent of formulation
    return build_initial_conditions(config, mesh)
end
