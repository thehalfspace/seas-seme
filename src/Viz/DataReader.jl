"""
    DataReader

Utilities for reading HDF5 simulation output data for visualization.
"""

using HDF5
using DelimitedFiles


"""
    read_timeseries_data(h5_file::String) -> (times, Vfmax)

Read global time series data from HDF5 output file.

# Arguments
- `h5_file::String`: Path to HDF5 output file

# Returns
- `times::Vector`: Simulation times [s]
- `Vfmax::Vector`: Maximum fault slip rate at each timestep [m/s]

# Example
```julia
times, Vfmax = read_timeseries_data("outputs/simulation.h5")
```
"""
function read_timeseries_data(h5_file::String)
    h5open(h5_file, "r") do file
        times = read(file, "timeseries/time")
        Vfmax = read(file, "timeseries/max_slip_rate")
        return times, Vfmax
    end
end


"""
    read_fault_geometry(h5_file::String) -> (x, y, depths_km)

Read fault geometry from HDF5 output file.

# Arguments
- `h5_file::String`: Path to HDF5 output file

# Returns
- `x::Vector`: Fault x coordinates [m]
- `y::Vector`: Fault y coordinates [m]
- `depths_km::Vector`: Fault depths [km]

# Example
```julia
x, y, depths_km = read_fault_geometry("outputs/simulation.h5")
```
"""
function read_fault_geometry(h5_file::String)
    h5open(h5_file, "r") do file
        x = read(file, "mesh/fault_x")
        y = read(file, "mesh/fault_y")
        depths_km = read(file, "mesh/fault_depths_km")
        return x, y, depths_km
    end
end


"""
    read_depth_timeseries(h5_file::String, depth::Float64) -> (time, slip, slip_rate, stress, state)

Read time series data for a specific depth from HDF5 file.

# Arguments
- `h5_file::String`: Path to HDF5 output file
- `depth::Float64`: Depth in km (e.g., 10.0 for 10 km)

# Returns
- `time::Vector`: Simulation times [s]
- `slip::Vector`: Cumulative slip [m]
- `slip_rate::Vector`: Slip rate [m/s]
- `stress::Vector`: Shear stress [MPa]
- `state::Vector`: State variable ψ

# Example
```julia
time, slip, slip_rate, stress, state = read_depth_timeseries("outputs/simulation.h5", 10.0)
```
"""
function read_depth_timeseries(h5_file::String, depth::Float64)
    depth_str = @sprintf("depth_%.1fkm", depth)

    h5open(h5_file, "r") do file
        time = read(file, "timeseries/time")
        slip = read(file, "timeseries/$depth_str/slip")
        slip_rate = read(file, "timeseries/$depth_str/slip_rate")
        stress = read(file, "timeseries/$depth_str/shear_stress")
        state = read(file, "timeseries/$depth_str/state_variable")

        return time, slip, slip_rate, stress, state
    end
end


"""
    read_initial_conditions(params_dir::String) -> (depths, a, b, Lc, σ_n, τ)

Read friction parameters and initial stresses from params directory.

# Arguments
- `params_dir::String`: Path to params directory (e.g., "data/simulation/params")

# Returns
- `depths::Vector`: Fault depths [km]
- `a::Vector`: Rate-state friction parameter a
- `b::Vector`: Rate-state friction parameter b
- `Lc::Vector`: Characteristic slip distance [m]
- `σ_n::Vector`: Normal stress [MPa]
- `τ::Vector`: Shear stress [MPa]

# Example
```julia
depths, a, b, Lc, σ_n, τ = read_initial_conditions("data/strike_slip_benchmark/params")
```
"""
function read_initial_conditions(params_dir::String)
    # Read friction parameters
    friction_file = joinpath(params_dir, "friction_parameters.dat")
    friction_data = readdlm(friction_file, skipstart=2)

    depths_friction = friction_data[:, 1]
    a = friction_data[:, 2]
    b = friction_data[:, 3]
    Lc = friction_data[:, 4]

    # Read initial stress
    stress_file = joinpath(params_dir, "initial_stress.dat")
    stress_data = readdlm(stress_file, skipstart=2)

    depths_stress = stress_data[:, 1]
    σ_n = stress_data[:, 2]
    τ = stress_data[:, 3]

    # Verify depths match
    if !all(isapprox.(depths_friction, depths_stress, atol=1e-6))
        @warn "Depth mismatch between friction and stress files"
    end

    return depths_friction, a, b, Lc, σ_n, τ
end


"""
    read_all_fault_timeseries(h5_file::String) -> (times, slip_rate_matrix, depths_km)

Read slip rate for all fault nodes to create heatmap data.

# Arguments
- `h5_file::String`: Path to HDF5 output file

# Returns
- `times::Vector`: Simulation times [s]
- `slip_rate_matrix::Matrix`: Slip rate matrix [n_depths × n_times] [m/s]
- `depths_km::Vector`: Fault depths [km]

# Example
```julia
times, slip_rates, depths = read_all_fault_timeseries("outputs/simulation.h5")
```
"""
function read_all_fault_timeseries(h5_file::String)
    h5open(h5_file, "r") do file
        # Read time
        times = read(file, "timeseries/time")

        # Read fault geometry
        depths_km = read(file, "mesh/fault_depths_km")
        n_depths = length(depths_km)
        n_times = length(times)

        # Get all depth groups
        ts_group = file["timeseries"]
        depth_groups = filter(k -> startswith(k, "depth_"), keys(ts_group))

        # Initialize matrix
        slip_rate_matrix = zeros(n_depths, n_times)

        # Read slip rate for each depth
        for depth_str in depth_groups
            # Get fault index from attributes
            depth_group = ts_group[depth_str]
            fault_idx = read(attributes(depth_group)["fault_index"])

            # Read slip rate
            slip_rate = read(depth_group["slip_rate"])

            # Store in matrix
            slip_rate_matrix[fault_idx, :] .= slip_rate
        end

        return times, slip_rate_matrix, depths_km
    end
end


"""
    get_available_depths(h5_file::String) -> Vector{Float64}

Get list of available depth locations in HDF5 file.

# Arguments
- `h5_file::String`: Path to HDF5 output file

# Returns
- `depths::Vector{Float64}`: Available depths [km]

# Example
```julia
depths = get_available_depths("outputs/simulation.h5")
```
"""
function get_available_depths(h5_file::String)
    h5open(h5_file, "r") do file
        ts_group = file["timeseries"]
        depth_groups = filter(k -> startswith(k, "depth_"), keys(ts_group))

        depths = Float64[]
        for depth_str in depth_groups
            depth_group = ts_group[depth_str]
            depth = read(attributes(depth_group)["depth_km"])
            push!(depths, depth)
        end

        return sort(depths)
    end
end


"""
    read_metadata(h5_file::String) -> Dict{String, Any}

Read simulation metadata from HDF5 file.

# Arguments
- `h5_file::String`: Path to HDF5 output file

# Returns
- `metadata::Dict{String, Any}`: Dictionary with simulation metadata

# Example
```julia
meta = read_metadata("outputs/simulation.h5")
println("Simulation: ", meta["simulation_name"])
```
"""
function read_metadata(h5_file::String)
    h5open(h5_file, "r") do file
        meta_group = file["metadata"]
        attrs = attributes(meta_group)

        metadata = Dict{String, Any}()
        for key in keys(attrs)
            metadata[key] = read(attrs[key])
        end

        return metadata
    end
end
