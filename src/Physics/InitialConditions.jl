"""
Initial condition setup for fault simulations.

Generates depth-dependent stress and friction parameters for rate-state friction.
"""

"""
    build_initial_conditions(config::PhysicsConfig, mesh::UnstructuredSEMesh) -> InitialConditions

Construct initial conditions for fault nodes from configuration and mesh.

# Arguments
- `config::PhysicsConfig`: Physics configuration
- `mesh::UnstructuredSEMesh`: Complete mesh with boundaries

# Returns
- `InitialConditions{Float64}`: Initial stresses and friction parameters

# Depth Zones
The depth-dependent friction structure (from surface to depth):
1. **0-2 km**: Velocity strengthening shallow zone
   - Transition from a-b = +0.003 to a-b = -0.0041
   - Prevents shallow seismicity
2. **2-12 km**: Velocity weakening seismogenic zone
   - a-b = -0.0041 (constant)
   - Earthquakes nucleate here
3. **12-17 km**: Transition zone
   - a-b transitions from -0.0041 to +0.015
4. **17-20 km**: Velocity strengthening deep zone
   - a-b = +0.0024 to +0.0047
   - Stable creep

# Stress Profile
- Normal stress: 10-50 MPa (increases with depth)
- Shear stress: 0.01-30 MPa (depth-dependent, slightly overstressed in seismogenic zone)
"""
function build_initial_conditions(
    config::PhysicsConfig,
    mesh::UnstructuredSEMesh{T}
) where {T<:AbstractFloat}
    # Extract fault coordinates (y-coordinates = depth along vertical fault)
    fault_y = mesh.boundaries.fault.coords[2, :]  # y-coordinates of fault nodes
    nfault = length(fault_y)

    @info "Generating initial conditions for $(nfault) fault nodes"

    # Initialize with default values
    a = 0.015 * ones(T, nfault)
    b = 0.019 * ones(T, nfault)
    σo = 50.0e6 * ones(T, nfault)
    τo = 22.0e6 * ones(T, nfault)
    Lc = 8.0e-3 * ones(T, nfault)

    # Depth zone boundaries (meters from surface)
    p1 = 2.0e3   # 2 km
    p2 = 12.0e3  # 12 km
    p3 = 17.0e3  # 17 km
    p4 = 20.0e3  # 20 km

    # Transform fault_y to depth coordinates (surface = 0)
    fault_depth = maximum(fault_y) - minimum(fault_y)
    if p4 > fault_depth
        @warn "Fault depth ($(fault_depth/1e3) km) is less than deepest zone boundary ($(p4/1e3) km)"
    end

    # Depth coordinate: top surface = 0 km
    depth = abs.(fault_y .- maximum(fault_y))

    # Identify nodes in each depth zone
    zone1 = findall(depth .<= p1)              # 0-2 km: shallow VS
    zone2 = findall(p1 .< depth .<= p2)        # 2-12 km: seismogenic VW
    zone3 = findall(p2 .< depth .<= p3)        # 12-17 km: transition
    zone4 = findall(p3 .< depth .<= p4)        # 17-20 km: deep VS
    zone5 = findall(depth .> p4)               # > 20 km: very deep VS

    # Set initial shear stress (τo) profile
    # Zone 1: Linear increase from 0.01 MPa to 30 MPa
    τo[zone1] .= 0.01e6 .+ 1e6 * (30.0 - 0.01) / (p1 - 0.0) .* (depth[zone1] .- 0.0)

    # Zone 2: Constant 30 MPa (slightly overstressed)
    τo[zone2] .= 30.0e6

    # Zone 3: Linear decrease from 30 MPa to 22.5 MPa
    τo[zone3] .= 30.0e6 .+ 1e6 * (22.5 - 30.0) / (p3 - p2) .* (depth[zone3] .- p2)

    # Zone 4: Constant 22.5 MPa
    τo[zone4] .= 22.5e6

    # Set initial normal stress (σo) profile
    # Zone 1: Linear increase from 10 MPa to 50 MPa
    σo[zone1] .= 10e6 .+ 1e6 * (50.0 - 10.0) / (p1 - 0.0) .* (depth[zone1] .- 0.0)

    # Zones 2-4: Constant 50 MPa
    σo[zone2] .= 50.0e6
    σo[zone3] .= 50.0e6
    σo[zone4] .= 50.0e6

    # Set rate-state friction parameters (a, b constant, vary a-b)
    # b = 0.019 (constant everywhere)

    if config.initial_conditions.velocity_strengthening_shallow
        # Zone 1: Velocity strengthening (a > b)
        # Transition from a-b = +0.003 to a-b = -0.0041
        a[zone1] .= b[zone1] .- 0.003 .+
                    (-0.0041 + 0.003) / (p1 - 0.0) .* (depth[zone1] .- 0.0)
    else
        # Alternative: Constant a-b = -0.0041 everywhere (no shallow VS)
        a[zone1] .= b[zone1] .- 0.0041
    end

    # Zone 2: Velocity weakening seismogenic zone (a < b)
    # a-b = -0.0041 (constant)
    a[zone2] .= b[zone2] .- 0.0041

    # Zone 3: Transition to velocity strengthening
    # a-b transitions from -0.0041 to +0.015
    a[zone3] .= b[zone3] .- 0.0041 .+
                (0.015 + 0.0041) / (p3 - p2) .* (depth[zone3] .- p2)

    # Zone 4: Deep velocity strengthening
    # a-b transitions from +0.015 to +0.0024
    a[zone4] .= b[zone4] .+ 0.015 .+
                (0.0024 - 0.0015) / (p4 - p3) .* (depth[zone4] .- p3)

    # Zone 5: Very deep velocity strengthening
    # a-b = +0.0047
    a[zone5] .= b[zone5] .+ 0.0047

    # Create friction law
    friction = RateStateFriction(
        a, b, Lc,
        config.reference_slip_rate,  # Vo
        config.reference_friction     # fo
    )

    # Create initial conditions
    ics = InitialConditions(σo, τo, friction)

    # Log summary statistics
    @info "Initial conditions summary:"
    @info "  a-b: min=$(minimum(a .- b)), max=$(maximum(a .- b))"
    @info "  Normal stress: $(minimum(σo)/1e6) - $(maximum(σo)/1e6) MPa"
    @info "  Shear stress: $(minimum(τo)/1e6) - $(maximum(τo)/1e6) MPa"
    @info "  Velocity strengthening shallow: $(config.initial_conditions.velocity_strengthening_shallow)"

    return ics
end

"""
    save_initial_conditions(ics::InitialConditions, fault_coords::Matrix, filepath::String)

Save initial conditions to text file for inspection.

# Arguments
- `ics::InitialConditions`: Initial conditions
- `fault_coords::Matrix`: Fault coordinates [2 x nfault]
- `filepath::String`: Output file path

# Output Format
Text file with columns:
- fault_depth (km)
- normal_stress (MPa)
- shear_stress (MPa)
- friction_a
- friction_b
- friction_Lc
"""
function save_initial_conditions(
    ics::InitialConditions,
    fault_coords::Matrix,
    filepath::String
)
    fault_y = fault_coords[2, :]
    depth_km = abs.(fault_y .- maximum(fault_y)) / 1e3  # km from surface

    # Create output matrix
    output = hcat(
        depth_km,
        ics.σn / 1e6,  # MPa
        ics.τs / 1e6,  # MPa
        ics.friction.a,
        ics.friction.b,
        ics.friction.Lc
    )

    # Write to file
    open(filepath, "w") do io
        write(io, "fault_depth(km) normal_stress(MPa) shear_stress(MPa) friction_a friction_b friction_Lc\n")
        writedlm(io, output)
    end

    @info "Initial conditions saved to: $filepath"
    return nothing
end
