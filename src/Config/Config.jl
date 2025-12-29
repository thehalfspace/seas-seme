"""
TOML configuration loading and validation.

This module provides functions to load simulation parameters from TOML files
and validate them before constructing typed configuration structs.
"""

"""
    load_config(toml_path::String) -> SimulationConfig

Load configuration from TOML file and return typed `SimulationConfig` struct.

# Arguments
- `toml_path::String`: Path to TOML configuration file

# Returns
- `SimulationConfig`: Validated configuration struct

# Throws
- `ErrorException`: If TOML file is invalid or parameters fail validation

# Example
```julia
config = load_config("examples/config/strike_slip_2d.toml")
```
"""
function load_config(toml_path::String)::SimulationConfig
    # Parse TOML file
    if !isfile(toml_path)
        error("Configuration file not found: $toml_path")
    end

    dict = TOML.parsefile(toml_path)

    # Get time conversion constants
    yr2sec = get(get(dict, "constants", Dict()), "yr2sec", 365.25 * 24 * 3600)
    day2sec = get(get(dict, "constants", Dict()), "day2sec", 24 * 3600)

    # Parse simulation section with time conversion
    sim_dict = dict["simulation"]
    total_time = if haskey(sim_dict, "total_time_years")
        sim_dict["total_time_years"] * yr2sec
    elseif haskey(sim_dict, "total_time")
        sim_dict["total_time"]
    else
        error("Must specify either total_time or total_time_years in [simulation]")
    end

    # Create simulation directory structure
    base_dir = sim_dict["output_dir"]
    sim_name = sim_dict["name"]
    sim_dir = joinpath(base_dir, sim_name)

    # Parse nested configurations (passing time constants for conversion)
    mesh_config = parse_mesh_config(dict["mesh"])
    physics_config = parse_physics_config(dict["physics"])
    solvers_config = parse_solvers_config(dict["solvers"], day2sec)
    output_config = parse_output_config(dict["output"], yr2sec)
    checkpoint_config = parse_checkpoint_config(dict["checkpointing"], sim_dir, yr2sec)

    # Construct top-level config
    sim_info = (
        name = sim_name,
        total_time = total_time,
        output_dir = sim_dir
    )

    config = SimulationConfig(
        sim_info,
        mesh_config,
        physics_config,
        solvers_config,
        output_config,
        checkpoint_config
    )

    # Validate entire configuration
    validate_config(config)

    return config
end

"""
    parse_mesh_config(dict::Dict) -> MeshConfig

Parse mesh configuration from TOML dict.
"""
function parse_mesh_config(dict::Dict)::MeshConfig
    return MeshConfig(
        dict["file"],
        dict["polynomial_degree"]
    )
end

"""
    parse_physics_config(dict::Dict) -> PhysicsConfig

Parse physics configuration from TOML dict.
"""
function parse_physics_config(dict::Dict)::PhysicsConfig
    ic_config = InitialConditionsConfig(
        get(dict["initial_conditions"], "type", "depth_dependent"),
        dict["initial_conditions"]["velocity_strengthening_shallow"],
        get(dict["initial_conditions"], "file", nothing)
    )

    return PhysicsConfig(
        dict["plate_velocity"],
        dict["reference_friction"],
        dict["reference_slip_rate"],
        dict["density"],
        dict["shear_velocity"],
        ic_config
    )
end

"""
    parse_solvers_config(dict::Dict, day2sec::Float64) -> SolverConfig

Parse solver configuration from TOML dict with time conversion.
"""
function parse_solvers_config(dict::Dict, day2sec::Float64)::SolverConfig
    qs_config = QuasistaticConfig(
        dict["quasistatic"]["tolerance"],
        dict["quasistatic"]["max_iterations"],
        dict["quasistatic"]["amg_max_levels"]
    )

    dyn_config = DynamicConfig(
        dict["dynamic"]["velocity_threshold_qs_to_dyn"],
        dict["dynamic"]["velocity_threshold_dyn_to_qs"]
    )

    # Convert dt_max from days to seconds if specified in days
    dt_max = if haskey(dict, "dt_max_days")
        dict["dt_max_days"] * day2sec
    elseif haskey(dict, "dt_max")
        dict["dt_max"]
    else
        error("Must specify either dt_max or dt_max_days in [solvers]")
    end

    return SolverConfig(
        dict["cfl_number"],
        dt_max,
        dict["dt_min_factor"],
        qs_config,
        dyn_config
    )
end

"""
    parse_output_config(dict::Dict, yr2sec::Float64) -> OutputConfig

Parse output configuration from TOML dict with time conversion.
"""
function parse_output_config(dict::Dict, yr2sec::Float64)::OutputConfig
    fields_config = OutputFieldsConfig(
        dict["fields"]["displacement"],
        dict["fields"]["velocity"],
        dict["fields"]["fault_slip"],
        dict["fields"]["fault_stress"],
        dict["fields"]["state_variable"]
    )

    # Convert snapshot interval from years to seconds if specified in years (backward compatibility)
    snapshot_interval = if haskey(dict, "snapshot_interval_years")
        dict["snapshot_interval_years"] * yr2sec
    elseif haskey(dict, "snapshot_interval")
        dict["snapshot_interval"]
    else
        1.0 * yr2sec  # Default to 1 year
    end

    # Parse snapshots configuration
    snapshots_config = if haskey(dict, "snapshots")
        snap_dict = dict["snapshots"]

        # Convert intervals
        qs_interval = if haskey(snap_dict, "quasistatic_interval_years")
            snap_dict["quasistatic_interval_years"] * yr2sec
        elseif haskey(snap_dict, "quasistatic_interval")
            snap_dict["quasistatic_interval"]
        else
            2.0 * yr2sec  # Default to 2 years
        end

        dyn_interval = get(snap_dict, "dynamic_interval", 0.1)  # Default to 0.1 seconds
        velocity_threshold = get(snap_dict, "velocity_threshold", 1e-3)  # Default to 1 mm/s

        SnapshotConfig(
            get(snap_dict, "enabled", true),
            qs_interval,
            dyn_interval,
            velocity_threshold
        )
    else
        # Default snapshot configuration
        SnapshotConfig(true, 2.0 * yr2sec, 0.1, 1e-3)
    end

    return OutputConfig(
        dict["format"],
        dict["timeseries_depths"],
        snapshot_interval,
        dict["log_interval"],
        fields_config,
        snapshots_config
    )
end

"""
    parse_checkpoint_config(dict::Dict, sim_dir::String, yr2sec::Float64) -> CheckpointConfig

Parse checkpointing configuration from TOML dict with time conversion.
"""
function parse_checkpoint_config(dict::Dict, sim_dir::String, yr2sec::Float64)::CheckpointConfig
    # Convert checkpoint interval from years to seconds if specified in years
    interval = if haskey(dict, "interval_years")
        dict["interval_years"] * yr2sec
    elseif haskey(dict, "interval")
        dict["interval"]
    else
        error("Must specify either interval or interval_years in [checkpointing]")
    end

    # Checkpoint directory inside simulation directory
    ckpt_dir = joinpath(sim_dir, "checkpoints")

    return CheckpointConfig(
        dict["enabled"],
        interval,
        dict["keep_last"],
        ckpt_dir
    )
end

"""
    validate_config(config::SimulationConfig)

Validate configuration parameters and throw errors if invalid.

# Validation checks:
- Positive time values
- Physical parameter ranges
- Solver parameters
- File paths existence
"""
function validate_config(config::SimulationConfig)
    # Simulation validation
    config.simulation.total_time > 0 || error("total_time must be positive")

    # Mesh validation
    config.mesh.polynomial_degree >= 1 || error("polynomial_degree must be >= 1")

    # Physics validation
    config.physics.plate_velocity > 0 || error("plate_velocity must be positive")
    config.physics.reference_friction > 0 || error("reference_friction must be positive")
    config.physics.reference_slip_rate > 0 || error("reference_slip_rate must be positive")
    config.physics.density > 0 || error("density must be positive")
    config.physics.shear_velocity > 0 || error("shear_velocity must be positive")

    # Solver validation
    0 < config.solvers.cfl_number <= 1 || error("cfl_number must be in (0, 1]")
    config.solvers.dt_max > 0 || error("dt_max must be positive")
    config.solvers.dt_min_factor > 0 || error("dt_min_factor must be positive")
    config.solvers.quasistatic.tolerance > 0 || error("QS tolerance must be positive")
    config.solvers.quasistatic.max_iterations > 0 || error("QS max_iterations must be positive")
    config.solvers.dynamic.velocity_threshold_qs_to_dyn > 0 ||
        error("velocity_threshold_qs_to_dyn must be positive")
    config.solvers.dynamic.velocity_threshold_dyn_to_qs > 0 ||
        error("velocity_threshold_dyn_to_qs must be positive")

    # Output validation
    config.output.format == "hdf5" || error("Only HDF5 output format currently supported")
    length(config.output.timeseries_depths) > 0 || error("Must specify at least one timeseries depth")
    config.output.snapshot_interval > 0 || error("snapshot_interval must be positive")
    config.output.log_interval > 0 || error("log_interval must be positive")

    # Snapshot validation
    if config.output.snapshots.enabled
        config.output.snapshots.quasistatic_interval > 0 ||
            error("quasistatic_interval must be positive")
        config.output.snapshots.dynamic_interval > 0 ||
            error("dynamic_interval must be positive")
        config.output.snapshots.velocity_threshold > 0 ||
            error("velocity_threshold must be positive")
    end

    # Checkpointing validation
    if config.checkpointing.enabled
        config.checkpointing.interval > 0 || error("checkpoint interval must be positive")
        config.checkpointing.keep_last > 0 || error("keep_last must be positive")
    end

    @info "Configuration validated successfully"
    return nothing
end
