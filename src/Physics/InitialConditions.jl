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
        a, b, Lc, σo,
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
    build_initial_conditions_from_csv(csv_path, mesh, physics_config, dip_angle) -> InitialConditions

Build initial conditions by reading friction/stress parameters from a CSV file
and interpolating onto mesh fault nodes.

The CSV file specifies parameters at discrete down-dip distances (x_d).
Shear stress τ⁰ and state variable ψ are computed from the steady-state
rate-state friction law (BP3-FD formulation, no radiation damping).

# Arguments
- `csv_path::String`: Path to CSV file with columns: x_d, a, b, Lc, sigma_n
- `mesh::UnstructuredSEMesh`: Complete mesh with fault boundary
- `physics_config::PhysicsConfig`: Physics configuration (for Vo, fo, Vpl)
- `dip_angle::Float64`: Fault dip angle in degrees

# CSV Format
```
x_d,a,b,Lc,sigma_n
0.0,0.010,0.015,0.008,50.0e6
15000.0,0.010,0.015,0.008,50.0e6
...
```
Where x_d is down-dip distance from free surface (meters).

# Notes
- Linear interpolation between CSV rows onto mesh fault nodes
- τ⁰ computed from BP3-FD eq. 22 (fully dynamic, no radiation damping):
  τ⁰ = σ_n * a_max * asinh(V_init/(2V₀) * exp((f₀ + b*ln(V₀/V_init))/a_max))
  where a_max = max(a) — spatially uniform τ⁰ seeds nucleation in VW zone
- V_init = plate_velocity (steady-state loading rate)
"""
function build_initial_conditions_from_csv(
    csv_path::String,
    mesh::UnstructuredSEMesh{T},
    physics_config::PhysicsConfig,
    dip_angle::Float64
) where {T<:AbstractFloat}
    @info "Loading initial conditions from CSV: $csv_path"

    # Read CSV file, skipping comment lines
    lines = filter(l -> !startswith(strip(l), '#') && !isempty(strip(l)),
                   readlines(csv_path))

    # Parse header and data
    header = strip.(split(lines[1], ','))
    expected = ["x_d", "a", "b", "Lc", "sigma_n"]
    header == expected || error("CSV header mismatch: expected $expected, got $header")

    nrows = length(lines) - 1
    csv_xd = zeros(T, nrows)
    csv_a  = zeros(T, nrows)
    csv_b  = zeros(T, nrows)
    csv_Lc = zeros(T, nrows)
    csv_σn = zeros(T, nrows)

    for (i, line) in enumerate(lines[2:end])
        vals = parse.(T, split(strip(line), ','))
        csv_xd[i] = vals[1]
        csv_a[i]  = vals[2]
        csv_b[i]  = vals[3]
        csv_Lc[i] = vals[4]
        csv_σn[i] = vals[5]
    end

    # Verify CSV is sorted by x_d
    issorted(csv_xd) || error("CSV x_d values must be sorted in ascending order")

    # Compute fault node down-dip distances from mesh coordinates
    fault_y = mesh.boundaries.fault.coords[2, :]
    nfault = length(fault_y)
    y_surface = maximum(fault_y)
    vertical_depth = abs.(fault_y .- y_surface)
    x_d_nodes = vertical_depth ./ T(sind(dip_angle))

    @info "  Fault nodes: $nfault"
    @info "  Down-dip distance range: $(minimum(x_d_nodes)/1e3) - $(maximum(x_d_nodes)/1e3) km"
    @info "  CSV x_d range: $(minimum(csv_xd)/1e3) - $(maximum(csv_xd)/1e3) km"

    # Linearly interpolate each parameter from CSV grid onto fault nodes
    a  = _interp_linear(csv_xd, csv_a,  x_d_nodes)
    b  = _interp_linear(csv_xd, csv_b,  x_d_nodes)
    Lc = _interp_linear(csv_xd, csv_Lc, x_d_nodes)
    σn = _interp_linear(csv_xd, csv_σn, x_d_nodes)

    # Compute τ⁰ from BP3-FD eq. 22 (fully dynamic steady-state)
    # τ⁰ = σ_n * a * asinh( V_init/(2*V₀) * exp((f₀ + b*ln(V₀/V_init)) / a) )
    # V_init = fault-tangential plate rate (after dip projection if dip-slip)
    Vo = T(physics_config.reference_slip_rate)
    fo = T(physics_config.reference_friction)
    V_plate = T(physics_config.plate_velocity)
    V_init = if physics_config.loading_direction == :dip_slip
        V_plate * T(cosd(dip_angle))
    else
        V_plate
    end

    # BP3 Eq. 22/25: τ⁰ uses a_max (not local a[i]) to create uniform stress
    # that is overstressed in VW zone, seeding nucleation
    a_max = maximum(a)
    τo = zeros(T, nfault)
    for i in 1:nfault
        arg = V_init / (2 * Vo) * exp((fo + b[i] * log(Vo / V_init)) / a_max)
        τo[i] = σn[i] * a_max * asinh(arg)
    end

    # Build friction and IC structs
    friction = RateStateFriction(a, b, Lc, σn, Vo, fo)
    ics = InitialConditions(σn, τo, friction)

    @info "Initial conditions from CSV:"
    @info "  V_init (plate rate on fault): $V_init m/s"
    @info "  a range: $(minimum(a)) - $(maximum(a))"
    @info "  b range: $(minimum(b)) - $(maximum(b))"
    @info "  a-b range: $(minimum(a .- b)) - $(maximum(a .- b))"
    @info "  σ_n range: $(minimum(σn)/1e6) - $(maximum(σn)/1e6) MPa"
    @info "  τ⁰ range: $(minimum(τo)/1e6) - $(maximum(τo)/1e6) MPa"

    return ics
end


"""
    _interp_linear(x_grid, y_grid, x_query) -> Vector

Piecewise linear interpolation. Clamps to boundary values for points
outside the grid range.
"""
function _interp_linear(x_grid::Vector{T}, y_grid::Vector{T}, x_query::Vector{T}) where T
    n = length(x_query)
    result = zeros(T, n)
    ng = length(x_grid)

    for i in 1:n
        xq = x_query[i]
        if xq <= x_grid[1]
            result[i] = y_grid[1]
        elseif xq >= x_grid[end]
            result[i] = y_grid[end]
        else
            # Find bracketing interval
            j = searchsortedlast(x_grid, xq)
            j = clamp(j, 1, ng - 1)
            t = (xq - x_grid[j]) / (x_grid[j+1] - x_grid[j])
            result[i] = y_grid[j] + t * (y_grid[j+1] - y_grid[j])
        end
    end

    return result
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
        ics.σo / 1e6,  # MPa
        ics.τo / 1e6,  # MPa
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
