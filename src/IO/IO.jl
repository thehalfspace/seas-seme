"""
IO module: Input/output operations for SEAS simulations.

This module provides:
- **HDF5Writer**: Streaming HDF5 output with compression
- **Checkpointing**: Save/load simulation state
- **Logging**: Structured event logging

# HDF5 Output Structure

```
simulation.h5
├── /metadata
│   ├── simulation_name
│   ├── creation_time
│   └── total_time
│
├── /mesh
│   ├── fault_x [nfault]
│   ├── fault_y [nfault]
│   └── fault_depths_km [nfault]
│
└── /timeseries
    ├── time [ntimesteps]
    ├── max_slip_rate [ntimesteps]
    │
    ├── /depth_0.0km
    │   ├── slip [ntimesteps]
    │   ├── slip_rate [ntimesteps]
    │   ├── shear_stress [ntimesteps]
    │   └── state_variable [ntimesteps]
    │
    └── /depth_20.0km
        └── ...
```

# Checkpointing

Checkpoints are saved as JLD2 files containing:
- Complete simulation state (u, v, a, fault variables)
- Solver objects (including AMG preconditioner)
- Configuration and metadata

# Usage Example

```julia
# Initialize I/O
io_manager = create_io_manager(config, mesh)
initialize_output!(io_manager, mesh, config, ics)

# Write time series
write_timestep!(io_manager, state, mesh, ics, params)

# Save checkpoint
save_checkpoint!(simulation, config)

# Load checkpoint
sim = load_checkpoint("checkpoint.jld2", mesh)

# Finalize
finalize_output!(io_manager)
```
"""

# Module exports
export HDF5OutputManager, create_io_manager
export initialize_output!, write_timestep!, write_snapshot!, should_write_snapshot, finalize_output!
export save_checkpoint!, load_checkpoint, find_latest_checkpoint
export setup_logging, close_logging, log_event
export setup_simulation_directories, save_initial_parameters
