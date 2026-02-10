"""
Configuration parameter types for SEAS-SEME simulations.

All parameters are loaded from TOML configuration files and validated.
"""

# ============================================================================
# Leaf types (no dependencies on other config types)
# ============================================================================

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

# ============================================================================
# Intermediate types (depend on leaf types)
# ============================================================================

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
- `formulation::Symbol`: Deformation formulation (:antiplane or :plane_strain)
- `dip_angle::Float64`: Fault dip angle in degrees (90 = vertical)
- `poisson_ratio::Float64`: Poisson's ratio ν (only used for plane_strain)
"""
struct PhysicsConfig
    plate_velocity::Float64
    reference_friction::Float64
    reference_slip_rate::Float64
    density::Float64
    shear_velocity::Float64
    initial_conditions::InitialConditionsConfig
    formulation::Symbol              # :antiplane or :plane_strain
    dip_angle::Float64               # Fault dip angle in degrees (90 = vertical)
    poisson_ratio::Float64           # Poisson's ratio ν (used for plane_strain)
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
    SnapshotConfig

Configuration for full fault snapshots at intervals.

# Fields
- `enabled::Bool`: Enable/disable snapshot output
- `quasistatic_interval::Float64`: Time interval during quasistatic periods (seconds)
- `dynamic_interval::Float64`: Time interval during dynamic events (seconds)
- `velocity_threshold::Float64`: Slip rate threshold to identify dynamic events (m/s)
"""
struct SnapshotConfig
    enabled::Bool
    quasistatic_interval::Float64
    dynamic_interval::Float64
    velocity_threshold::Float64
end

"""
    OutputConfig

Output configuration (HDF5, time series, snapshots).

# Fields
- `format::String`: Output format ("hdf5")
- `timeseries_depths::Vector{Float64}`: Depths (km) for time series output
- `snapshot_interval::Float64`: Time interval for full-field snapshots (seconds) [DEPRECATED - use snapshots config]
- `log_interval::Int`: Iteration interval for console logging
- `fields::OutputFieldsConfig`: Which fields to output
- `snapshots::SnapshotConfig`: Configuration for full fault snapshots
"""
struct OutputConfig
    format::String
    timeseries_depths::Vector{Float64}
    snapshot_interval::Float64  # Keep for backward compatibility
    log_interval::Int
    fields::OutputFieldsConfig
    snapshots::SnapshotConfig
end

# ============================================================================
# Top-level type (depends on all other types)
# ============================================================================

"""
    SimulationConfig

Top-level configuration struct containing all simulation parameters.

# Fields
- `simulation::NamedTuple`: Simulation metadata (name, total_time, output_dir)
- `mesh::MeshConfig`: Mesh configuration
- `physics::PhysicsConfig`: Physics parameters
- `solvers::SolverConfig`: Solver settings
- `output::OutputConfig`: Output configuration
- `checkpointing::CheckpointConfig`: Checkpointing settings
"""
struct SimulationConfig
    simulation::NamedTuple{(:name, :total_time, :output_dir), Tuple{String, Float64, String}}
    mesh::MeshConfig
    physics::PhysicsConfig
    solvers::SolverConfig
    output::OutputConfig
    checkpointing::CheckpointConfig
end
