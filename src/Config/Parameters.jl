"""
Configuration parameter types for SEAS-SEME simulations.

All parameters are loaded from TOML configuration files and validated.
"""

"""
    SimulationConfig

Top-level configuration struct containing all simulation parameters.

# Fields
- `name::String`: Simulation name (used for output directories)
- `total_time::Float64`: Total simulation time (seconds)
- `output_dir::String`: Base directory for output files
- `mesh::MeshConfig`: Mesh configuration
- `physics::PhysicsConfig`: Physics parameters
- `solvers::SolverConfig`: Solver settings
- `output::OutputConfig`: Output configuration
- `checkpointing::CheckpointConfig`: Checkpointing settings
"""
struct SimulationConfig
    name::String
    total_time::Float64
    output_dir::String
    mesh::MeshConfig
    physics::PhysicsConfig
    solvers::SolverConfig
    output::OutputConfig
    checkpointing::CheckpointConfig
end

"""
    MeshConfig

Mesh generation and discretization parameters.

# Fields
- `file::String`: Path to mesh file (.mesh format from HOHQMesh)
- `polynomial_degree::Int`: Polynomial degree for spectral elements (p)
"""
struct MeshConfig
    file::String
    polynomial_degree::Int
end

"""
    PhysicsConfig

Physical parameters for the simulation.

# Fields
- `plate_velocity::Float64`: Plate loading rate (m/s)
- `reference_friction::Float64`: Reference friction coefficient (dimensionless)
- `reference_slip_rate::Float64`: Reference slip rate Vo (m/s)
- `density::Float64`: Material density (kg/m³)
- `shear_velocity::Float64`: Shear wave velocity (m/s)
- `initial_conditions::InitialConditionsConfig`: Initial condition parameters
"""
struct PhysicsConfig
    plate_velocity::Float64
    reference_friction::Float64
    reference_slip_rate::Float64
    density::Float64
    shear_velocity::Float64
    initial_conditions::InitialConditionsConfig
end

"""
    InitialConditionsConfig

Configuration for initial conditions (stress, friction parameters).

# Fields
- `type::String`: Type of initial conditions ("depth_dependent" or "custom")
- `velocity_strengthening_shallow::Bool`: Use velocity-strengthening friction in shallow zone
- `file::Union{String,Nothing}`: Optional path to custom IC file
"""
struct InitialConditionsConfig
    type::String
    velocity_strengthening_shallow::Bool
    file::Union{String,Nothing}
end

"""
    SolverConfig

Solver configuration (quasistatic, dynamic, timestepping).

# Fields
- `cfl_number::Float64`: CFL number for stability (typically 0.5-0.6)
- `dt_max::Float64`: Maximum timestep (seconds)
- `dt_min_factor::Float64`: Factor for minimum timestep (dt_min = factor * dx/vs)
- `quasistatic::QuasistaticConfig`: Quasistatic solver settings
- `dynamic::DynamicConfig`: Dynamic solver settings
"""
struct SolverConfig
    cfl_number::Float64
    dt_max::Float64
    dt_min_factor::Float64
    quasistatic::QuasistaticConfig
    dynamic::DynamicConfig
end

"""
    QuasistaticConfig

Quasistatic solver parameters (AMG-preconditioned CG).

# Fields
- `tolerance::Float64`: Convergence tolerance for linear solver
- `max_iterations::Int`: Maximum CG iterations
- `amg_max_levels::Int`: Maximum AMG levels for preconditioner
"""
struct QuasistaticConfig
    tolerance::Float64
    max_iterations::Int
    amg_max_levels::Int
end

"""
    DynamicConfig

Dynamic solver parameters (explicit time integration).

# Fields
- `velocity_threshold_qs_to_dyn::Float64`: Slip rate threshold to switch QS → dynamic (m/s)
- `velocity_threshold_dyn_to_qs::Float64`: Slip rate threshold to switch dynamic → QS (m/s)
"""
struct DynamicConfig
    velocity_threshold_qs_to_dyn::Float64
    velocity_threshold_dyn_to_qs::Float64
end

"""
    OutputConfig

Output configuration (HDF5, time series, snapshots).

# Fields
- `format::String`: Output format ("hdf5")
- `timeseries_depths::Vector{Float64}`: Depths (km) for time series output
- `snapshot_interval::Float64`: Time interval for full-field snapshots (seconds)
- `log_interval::Int`: Iteration interval for console logging
- `fields::OutputFieldsConfig`: Which fields to output
"""
struct OutputConfig
    format::String
    timeseries_depths::Vector{Float64}
    snapshot_interval::Float64
    log_interval::Int
    fields::OutputFieldsConfig
end

"""
    OutputFieldsConfig

Flags for which fields to include in output.

# Fields
- `displacement::Bool`: Output displacement field
- `velocity::Bool`: Output velocity field
- `fault_slip::Bool`: Output fault slip
- `fault_stress::Bool`: Output fault shear stress
- `state_variable::Bool`: Output rate-state variable (ψ or θ)
"""
struct OutputFieldsConfig
    displacement::Bool
    velocity::Bool
    fault_slip::Bool
    fault_stress::Bool
    state_variable::Bool
end

"""
    CheckpointConfig

Checkpointing configuration.

# Fields
- `enabled::Bool`: Enable/disable checkpointing
- `interval::Float64`: Time interval between checkpoints (seconds)
- `keep_last::Int`: Number of checkpoints to keep (oldest deleted)
- `directory::String`: Directory for checkpoint files
"""
struct CheckpointConfig
    enabled::Bool
    interval::Float64
    keep_last::Int
    directory::String
end
