"""
    HDF5Writer

Streaming HDF5 output for earthquake cycle simulations.

Features:
- Extensible datasets for time series
- GZIP compression (level 4)
- Chunked storage for efficient access
- Metadata storage in attributes
- Zero-copy writes during simulation
"""

using HDF5
using Dates


"""
    HDF5OutputManager

Manager for HDF5 output during simulation.

# Fields
- `file::HDF5.File`: Open HDF5 file handle
- `output_indices::Vector{Int}`: Fault node indices for output
- `depths::Vector{Float64}`: Fault depths for output locations [km]
- `current_index::Int`: Current write index in datasets
- `write_interval::Int`: Write every N iterations
"""
mutable struct HDF5OutputManager
    file::HDF5.File
    output_indices::Vector{Int}
    depths::Vector{Float64}
    current_index::Int
    write_interval::Int
end


"""
    setup_simulation_directories(config::SimulationConfig, config_path::String)

Create simulation directory structure and copy config file.

# Creates:
- `simulation_name/`
  - `config.toml` (copy of input config)
  - `params/` (for friction parameters, stresses)
  - `outputs/` (for HDF5 files)
  - `checkpoints/` (for checkpoint files)
"""
function setup_simulation_directories(config::SimulationConfig, config_path::String)
    sim_dir = config.simulation.output_dir

    # Create directory structure
    mkpath(joinpath(sim_dir, "params"))
    mkpath(joinpath(sim_dir, "outputs"))
    mkpath(joinpath(sim_dir, "checkpoints"))

    # Copy config file to simulation directory
    config_dest = joinpath(sim_dir, "config.toml")
    if isfile(config_path)
        cp(config_path, config_dest, force=true)
        @info "Configuration copied to simulation directory" path=config_dest
    end

    @info "Simulation directories created" directory=sim_dir

    return nothing
end


"""
    create_io_manager(config::SimulationConfig, mesh) -> HDF5OutputManager

Create HDF5 output manager and initialize output file.

# Arguments
- `config::SimulationConfig`: Simulation configuration
- `mesh`: Spectral element mesh

# Returns
- `HDF5OutputManager`: Manager ready for time series writing
"""
function create_io_manager(config::SimulationConfig, mesh)
    # Output directory is simulation_name/outputs/
    output_dir = joinpath(config.simulation.output_dir, "outputs")
    mkpath(output_dir)

    # Output file path
    filepath = joinpath(output_dir, "$(config.simulation.name).h5")

    # Remove existing file if present
    isfile(filepath) && rm(filepath)

    # Create HDF5 file
    file = h5open(filepath, "w")

    # Find output locations (indices matching requested depths)
    fault_coords = mesh.boundaries.fault.coords
    fault_depths_km = -fault_coords[2, :] / 1000  # Convert to km, positive down

    output_indices = Int[]
    output_depths = Float64[]

    for target_depth in config.output.timeseries_depths
        # Find closest fault node to target depth
        idx = argmin(abs.(fault_depths_km .- target_depth))
        push!(output_indices, idx)
        push!(output_depths, fault_depths_km[idx])
    end

    @info "HDF5 output initialized" filepath=filepath n_locations=length(output_indices)

    return HDF5OutputManager(file, output_indices, output_depths, 0, 1)
end


"""
    save_initial_parameters(config::SimulationConfig, ics, mesh)

Save friction parameters and initial stresses to params/ directory.

# Saves:
- `params/friction_parameters.dat`: a, b, Lc along fault
- `params/initial_stress.dat`: σ_n, τ along fault
- `params/fault_coordinates.dat`: x, y coordinates
"""
function save_initial_parameters(config::SimulationConfig, ics, mesh)
    params_dir = joinpath(config.simulation.output_dir, "params")
    mkpath(params_dir)

    fault_coords = mesh.boundaries.fault.coords

    # Save friction parameters (a, b, Lc)
    friction_file = joinpath(params_dir, "friction_parameters.dat")
    open(friction_file, "w") do io
        println(io, "# Friction parameters along fault")
        println(io, "# Depth(km)  a  b  Lc(m)")
        for i in eachindex(ics.friction.a)
            depth_km = -fault_coords[2, i] / 1000
            @printf(io, "%.6f  %.6e  %.6e  %.6e\n",
                   depth_km, ics.friction.a[i], ics.friction.b[i], ics.friction.Lc[i])
        end
    end

    # Save initial stress (σ_n, τ)
    stress_file = joinpath(params_dir, "initial_stress.dat")
    open(stress_file, "w") do io
        println(io, "# Initial stress along fault")
        println(io, "# Depth(km)  σ_n(MPa)  τ(MPa)")
        for i in eachindex(ics.σo)
            depth_km = -fault_coords[2, i] / 1000
            @printf(io, "%.6f  %.6e  %.6e\n",
                   depth_km, ics.σo[i]/1e6, ics.τo[i]/1e6)
        end
    end

    # Save fault coordinates
    coords_file = joinpath(params_dir, "fault_coordinates.dat")
    open(coords_file, "w") do io
        println(io, "# Fault coordinates")
        println(io, "# x(m)  y(m)  Depth(km)")
        for i in 1:size(fault_coords, 2)
            depth_km = -fault_coords[2, i] / 1000
            @printf(io, "%.6e  %.6e  %.6f\n",
                   fault_coords[1, i], fault_coords[2, i], depth_km)
        end
    end

    @info "Initial parameters saved" directory=params_dir

    return nothing
end


"""
    initialize_output!(io::HDF5OutputManager, mesh, config, ics)

Initialize HDF5 file structure with metadata and extensible datasets.

# Creates:
- `/metadata` group with simulation metadata
- `/mesh` group with mesh information
- `/timeseries` group with extensible datasets for each output location
"""
function initialize_output!(io::HDF5OutputManager, mesh, config, ics)
    file = io.file

    # === Metadata ===
    meta = create_group(file, "metadata")
    attributes(meta)["simulation_name"] = config.simulation.name
    attributes(meta)["creation_time"] = string(now())
    attributes(meta)["total_time"] = config.simulation.total_time
    attributes(meta)["julia_version"] = string(VERSION)
    # TODO: Add git commit if available

    # === Mesh ===
    mesh_group = create_group(file, "mesh")
    fault_coords = mesh.boundaries.fault.coords
    mesh_group["fault_x"] = fault_coords[1, :]
    mesh_group["fault_y"] = fault_coords[2, :]
    mesh_group["fault_depths_km"] = -fault_coords[2, :] / 1000

    # === Time series datasets ===
    ts_group = create_group(file, "timeseries")

    # Global maximum slip rate
    create_extensible_dataset(ts_group, "time", Float64)
    create_extensible_dataset(ts_group, "max_slip_rate", Float64)

    # Per-location time series
    for (i, depth) in enumerate(io.depths)
        depth_str = @sprintf("depth_%.1fkm", depth)

        # Check if group already exists; if so, skip (can happen with duplicate depths)
        if haskey(ts_group, depth_str)
            @warn "Skipping duplicate depth group" depth_str depth
            continue
        end

        depth_group = create_group(ts_group, depth_str)

        # Store location metadata
        attributes(depth_group)["depth_km"] = depth
        attributes(depth_group)["fault_index"] = io.output_indices[i]

        # Create extensible datasets
        create_extensible_dataset(depth_group, "slip", Float64)
        create_extensible_dataset(depth_group, "slip_rate", Float64)
        create_extensible_dataset(depth_group, "shear_stress", Float64)
        create_extensible_dataset(depth_group, "state_variable", Float64)
    end

    # Flush to disk
    flush(file)

    @info "HDF5 structure initialized" n_depths=length(io.depths)

    return nothing
end


"""
    create_extensible_dataset(parent, name, dtype; chunk_size=1000)

Create extensible 1D dataset with compression.

# Arguments
- `parent`: HDF5 group
- `name::String`: Dataset name
- `dtype::Type`: Data type
- `chunk_size::Int`: Chunk size for storage (default: 1000)

# Returns
- `HDF5.Dataset`: Extensible dataset
"""
function create_extensible_dataset(parent, name::String, dtype::Type;
                                  chunk_size::Int=1000)
    # Create extensible dataset with GZIP compression
    dset = create_dataset(parent, name, dtype,
                         ((0,), (-1,)),  # Initial size 0, unlimited max
                         chunk=(chunk_size,),
                         deflate=4)  # GZIP compression level 4

    return dset
end


"""
    write_timestep!(io::HDF5OutputManager, state, mesh, ics, params)

Write current timestep data to HDF5 file.

# Arguments
- `io::HDF5OutputManager`: Output manager
- `state::SimulationState`: Current state
- `mesh`: Mesh with boundary info
- `ics`: Initial conditions
- `params`: Simulation parameters

# Notes
Extends datasets and appends new data. Automatically flushes every 100 writes.
"""
function write_timestep!(io::HDF5OutputManager, state, mesh, ics, params)
    file = io.file
    io.current_index += 1
    idx = io.current_index

    fault_id = mesh.boundaries.fault.node_ids

    # === Global time series ===
    extend_and_write!(file["timeseries/time"], idx, state.time)

    Vf_max = maximum(abs.(state.Vf))
    extend_and_write!(file["timeseries/max_slip_rate"], idx, Vf_max)

    # === Per-location time series ===
    for (i, fault_idx) in enumerate(io.output_indices)
        depth_str = @sprintf("depth_%.1fkm", io.depths[i])
        depth_path = "timeseries/$depth_str"

        # Compute quantities
        slip = 2 * state.u[fault_id[fault_idx]] + params.Vpl * state.time
        slip_rate = state.Vf[fault_idx]
        shear_stress = (state.τf[fault_idx] + ics.τo[fault_idx]) / 1e6  # MPa
        psi = state.ψ[fault_idx]

        # Write to datasets
        extend_and_write!(file["$depth_path/slip"], idx, slip)
        extend_and_write!(file["$depth_path/slip_rate"], idx, slip_rate)
        extend_and_write!(file["$depth_path/shear_stress"], idx, shear_stress)
        extend_and_write!(file["$depth_path/state_variable"], idx, psi)
    end

    # Periodic flush to disk (every 100 writes)
    if mod(idx, 100) == 0
        flush(file)
    end

    return nothing
end


"""
    extend_and_write!(dset::HDF5.Dataset, index::Int, value)

Extend dataset and write value at given index.

# Arguments
- `dset::HDF5.Dataset`: Dataset to write to
- `index::Int`: Index to write (1-based)
- `value`: Value to write
"""
function extend_and_write!(dset::HDF5.Dataset, index::Int, value)
    # Extend dataset to accommodate new index
    HDF5.set_extent_dims(dset, (index,))

    # Write value
    dset[index] = value

    return nothing
end


"""
    finalize_output!(io::HDF5OutputManager)

Finalize output and close HDF5 file.
"""
function finalize_output!(io::HDF5OutputManager)
    flush(io.file)
    close(io.file)

    @info "HDF5 output finalized" total_writes=io.current_index

    return nothing
end
