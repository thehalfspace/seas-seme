# Summary of Recent Changes

## 1. Improved Output Directory Structure

Each simulation now creates an organized directory structure:

```
data/
  simulation_name/
    config.toml                    # Copy of configuration used
    params/
      friction_parameters.dat      # a, b, Lc along fault
      initial_stress.dat           # σ_n, τ along fault
      fault_coordinates.dat        # x, y coordinates
    outputs/
      simulation_name.h5           # HDF5 time series data
    checkpoints/
      checkpoint_iter_*.jld2       # Checkpoint files
```

### Benefits:
- **Self-contained**: Each simulation has all its files in one place
- **Reproducible**: Config file is copied, so you know exactly what was run
- **Organized**: Clear separation of parameters, outputs, and checkpoints
- **Easy to compare**: Different simulations in different folders

## 2. Human-Readable Time Units in Config

### Old Format (seconds everywhere):
```toml
[simulation]
total_time = 6.31152e8  # Hard to read - what is this in years?

[solvers]
dt_max = 8.64e6  # Hard to read - what is this in days?

[checkpointing]
interval = 3.1557e7  # Hard to read - what is this?
```

### New Format (human-readable):
```toml
[constants]
yr2sec = 31557600.0   # Defined once
day2sec = 86400.0     # Defined once

[simulation]
total_time_years = 20.0  # Clear: 20 years!

[solvers]
dt_max_days = 100.0  # Clear: 100 days!

[checkpointing]
interval_years = 1.0  # Clear: 1 year!
```

### Backward Compatibility:
The parser still accepts the old format (`total_time`, `dt_max`, `interval` in seconds).
You can mix and match:
- Use `total_time` OR `total_time_years`
- Use `dt_max` OR `dt_max_days`
- Use `interval` OR `interval_years`

## 3. Automatic Parameter Saving

When a simulation starts, it now automatically saves:

### `params/friction_parameters.dat`:
```
# Friction parameters along fault
# Depth(km)  a  b  Lc(m)
0.000000  1.500000e-02  1.900000e-02  8.000000e-03
2.500000  1.200000e-02  1.600000e-02  8.000000e-03
...
```

### `params/initial_stress.dat`:
```
# Initial stress along fault
# Depth(km)  σ_n(MPa)  τ(MPa)
0.000000  2.670000e+01  1.602000e+01
2.500000  6.675000e+01  4.005000e+01
...
```

### `params/fault_coordinates.dat`:
```
# Fault coordinates
# x(m)  y(m)  Depth(km)
4.000000e+04  0.000000e+00  0.000000
4.000000e+04  -2.500000e+03  2.500000
...
```

## 4. Config File Automatically Copied

The exact config file used is copied to `simulation_name/config.toml`, ensuring:
- **Reproducibility**: You know exactly what parameters were used
- **Documentation**: No guessing what settings were used months later
- **Version control**: Each simulation has its own config snapshot

## Usage Examples

### Running a Simulation:
```bash
julia --project=. examples/scripts/run_simulation.jl examples/config/strike_slip_2d.toml
```

This creates:
```
data/strike_slip_benchmark/
├── config.toml
├── params/
│   ├── friction_parameters.dat
│   ├── initial_stress.dat
│   └── fault_coordinates.dat
├── outputs/
│   └── strike_slip_benchmark.h5
└── checkpoints/
    └── (checkpoints created during run)
```

### Examining Results:
```bash
# Check what parameters were used
cat data/strike_slip_benchmark/config.toml

# Plot friction parameters
python plot_friction.py data/strike_slip_benchmark/params/friction_parameters.dat

# Analyze HDF5 output
python analyze.py data/strike_slip_benchmark/outputs/strike_slip_benchmark.h5
```

### Restarting:
```bash
# Restart from latest checkpoint
julia --project=. examples/scripts/restart.jl --latest data/strike_slip_benchmark/checkpoints/
```

## Migration Guide

### Updating Existing Configs:

1. **Add constants section** (optional but recommended):
```toml
[constants]
yr2sec = 31557600.0
day2sec = 86400.0
```

2. **Convert time values** (optional but recommended):
```toml
# Old:
total_time = 6.31152e8

# New:
total_time_years = 20.0
```

3. **Remove checkpoint directory** (now auto-generated):
```toml
# Old:
[checkpointing]
directory = "./checkpoints"

# New: (directory field removed, auto-generated)
[checkpointing]
interval_years = 1.0
keep_last = 5
```

### Notes:
- Old config files still work! Backward compatible.
- Just add the new fields when you want human-readable times.
- The simulation directory structure is created automatically.

## What's Next

These changes make it ready for:
1. **Running multiple simulations**: Each gets its own directory
2. **Parameter sweeps**: Easy to organize many runs
3. **Reproducibility**: Everything needed to reproduce a run is in one place
4. **Collaboration**: Easy to share simulation results with all metadata
