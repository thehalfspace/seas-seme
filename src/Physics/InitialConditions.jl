"""
Initial condition setup for fault simulations.

Generates stress and friction parameters from a CSV file for rate-state friction.
"""

"""
    build_initial_conditions(csv_path, mesh, physics_config, dip_angle) -> InitialConditions

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
function build_initial_conditions(
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
