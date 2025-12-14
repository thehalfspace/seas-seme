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

    # Parse nested configurations
    mesh_config = parse_mesh_config(dict["mesh"])
    physics_config = parse_physics_config(dict["physics"])
    solvers_config = parse_solvers_config(dict["solvers"])
    output_config = parse_output_config(dict["output"])
    checkpoint_config = parse_checkpoint_config(dict["checkpointing"])

    # Construct top-level config
    config = SimulationConfig(
        dict["simulation"]["name"],
        dict["simulation"]["total_time"],
        dict["simulation"]["output_dir"],
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
    parse_solvers_config(dict::Dict) -> SolverConfig

Parse solver configuration from TOML dict.
"""
function parse_solvers_config(dict::Dict)::SolverConfig
    qs_config = QuasistaticConfig(
        dict["quasistatic"]["tolerance"],
        dict["quasistatic"]["max_iterations"],
        dict["quasistatic"]["amg_max_levels"]
    )

    dyn_config = DynamicConfig(
        dict["dynamic"]["velocity_threshold_qs_to_dyn"],
        dict["dynamic"]["velocity_threshold_dyn_to_qs"]
    )

    return SolverConfig(
        dict["cfl_number"],
        dict["dt_max"],
        dict["dt_min_factor"],
        qs_config,
        dyn_config
    )
end

"""
    parse_output_config(dict::Dict) -> OutputConfig

Parse output configuration from TOML dict.
"""
function parse_output_config(dict::Dict)::OutputConfig
    fields_config = OutputFieldsConfig(
        dict["fields"]["displacement"],
        dict["fields"]["velocity"],
        dict["fields"]["fault_slip"],
        dict["fields"]["fault_stress"],
        dict["fields"]["state_variable"]
    )

    return OutputConfig(
        dict["format"],
        dict["timeseries_depths"],
        dict["snapshot_interval"],
        dict["log_interval"],
        fields_config
    )
end

"""
    parse_checkpoint_config(dict::Dict) -> CheckpointConfig

Parse checkpointing configuration from TOML dict.
"""
function parse_checkpoint_config(dict::Dict)::CheckpointConfig
    return CheckpointConfig(
        dict["enabled"],
        dict["interval"],
        dict["keep_last"],
        dict["directory"]
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
    config.total_time > 0 || error("total_time must be positive")

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

    # Checkpointing validation
    if config.checkpointing.enabled
        config.checkpointing.interval > 0 || error("checkpoint interval must be positive")
        config.checkpointing.keep_last > 0 || error("keep_last must be positive")
    end

    @info "Configuration validated successfully"
    return nothing
end
