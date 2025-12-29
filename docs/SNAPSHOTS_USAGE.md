# Full Fault Snapshots and Slip Contour Plotting

This document explains how to use the snapshot functionality to create cumulative slip contour plots.

## Overview

The snapshot feature saves full fault profiles (slip, slip rate, stress, state variable) at regular intervals during the simulation. This allows you to create space-time contour plots showing the evolution of fault quantities over the entire earthquake cycle.

## Configuration

Add the `[output.snapshots]` section to your TOML configuration file:

```toml
[output.snapshots]
enabled = true                      # Enable full fault snapshots
quasistatic_interval_years = 2.0    # Snapshot interval during interseismic (years)
dynamic_interval = 0.1              # Snapshot interval during earthquakes (seconds)
velocity_threshold = 1.0e-3         # Threshold to identify dynamic events (1 mm/s)
```

### Parameters

- **`enabled`**: Turn snapshots on/off
- **`quasistatic_interval_years`**: How often to save snapshots during slow (quasi-static) periods, in years
- **`dynamic_interval`**: How often to save snapshots during fast (dynamic/seismic) events, in seconds
- **`velocity_threshold`**: Maximum slip rate threshold (m/s) to distinguish between quasi-static and dynamic regimes

## Output Format

Snapshots are saved in the HDF5 output file under `/snapshots`:

```
simulation.h5
└── /snapshots
    ├── times [n_snapshots]                      # Snapshot times (s)
    ├── slip [n_fault_nodes × n_snapshots]       # Cumulative slip (m)
    ├── slip_rate [n_fault_nodes × n_snapshots]  # Slip rate (m/s)
    ├── shear_stress [n_fault_nodes × n_snapshots] # Shear stress (MPa)
    └── state_variable [n_fault_nodes × n_snapshots] # State variable ψ
```

## Reading Snapshot Data

### Read all snapshots

```julia
using SEAS_SEME.Viz

# Read all snapshot data
times, slip, slip_rate, stress, state, depths_km = read_snapshots("data/strike_slip_2d/outputs/strike_slip_2d.h5")

# times: Vector of snapshot times (seconds)
# slip: Matrix [n_fault_nodes × n_snapshots] cumulative slip (m)
# slip_rate: Matrix [n_fault_nodes × n_snapshots] slip rate (m/s)
# stress: Matrix [n_fault_nodes × n_snapshots] shear stress (MPa)
# state: Matrix [n_fault_nodes × n_snapshots] state variable
# depths_km: Vector of fault node depths (km)
```

### Read single snapshot

```julia
# Read snapshot at index 10
time, slip, slip_rate, stress, state = read_snapshot_at_index("data/strike_slip_2d/outputs/strike_slip_2d.h5", 10)
```

### Get snapshot configuration

```julia
config = get_snapshot_config("data/strike_slip_2d/outputs/strike_slip_2d.h5")
println("Quasi-static interval: ", config["quasistatic_interval"])
```

## Creating Slip Contour Plots

### Quick plot - all contours

```julia
using SEAS_SEME.Viz

# Generate all contour plots in one call
plot_all_contours("data/strike_slip_2d/outputs/strike_slip_2d.h5",
                 output_dir="data/strike_slip_2d/plots",
                 time_units="years",
                 max_depth=20.0)  # Optional: limit depth
```

This creates three plots:
- `slip_contours.png` - Cumulative slip vs depth and time
- `slip_rate_contours.png` - Slip rate vs depth and time (log scale)
- `stress_contours.png` - Shear stress vs depth and time

### Individual contour plots

#### Cumulative slip contours

```julia
plot_slip_contours("data/strike_slip_2d/outputs/strike_slip_2d.h5",
                  output_file="slip_contours.png",
                  time_units="years",      # "years" or "seconds"
                  max_depth=20.0,          # km (optional)
                  colormap=:viridis,       # color scheme
                  levels=20)               # number of contour levels
```

#### Slip rate contours

```julia
plot_slip_rate_contours("data/strike_slip_2d/outputs/strike_slip_2d.h5",
                       output_file="slip_rate_contours.png",
                       time_units="years",
                       max_depth=20.0,
                       colormap=:hot,
                       log_scale=true)      # logarithmic color scale
```

#### Shear stress contours

```julia
plot_stress_contours("data/strike_slip_2d/outputs/strike_slip_2d.h5",
                    output_file="stress_contours.png",
                    time_units="years",
                    max_depth=20.0,
                    colormap=:balance)
```

## Complete Workflow Example

```julia
using SEAS_SEME

# 1. Load configuration with snapshots enabled
config = load_config("config/strike_slip_2d.toml")

# 2. Run simulation (snapshots will be saved automatically)
sim = setup_simulation(config, "config/strike_slip_2d.toml")
run!(sim)

# 3. Generate all plots
using SEAS_SEME.Viz

plot_all_contours("data/strike_slip_2d/outputs/strike_slip_2d.h5",
                 output_dir="data/strike_slip_2d/plots",
                 time_units="years")
```

## Comparison with Spear Format

The snapshot functionality is similar to Spear's output format:

| Spear | SEAS-SEME | Description |
|-------|-----------|-------------|
| `delfyr.out` | Snapshots with `quasistatic_interval` | Slip profiles every N years (interseismic) |
| `delfsec.out` | Snapshots with `dynamic_interval` | Slip profiles every N seconds (coseismic) |
| Text file | HDF5 `/snapshots` group | Storage format |

**Advantages of SEAS-SEME snapshots:**
- Single HDF5 file (not separate text files)
- Compressed storage (GZIP level 4)
- Automatic adaptive sampling (dense during events, sparse during quasi-static)
- Stores slip rate, stress, and state variable in addition to slip
- Easy to read in Julia, Python, or MATLAB

## Performance Considerations

- **Storage**: Each snapshot stores 4 arrays of size `n_fault_nodes`
- **Recommendation**: Use longer intervals for long simulations (e.g., 2-5 years for interseismic, 0.1-0.5s for dynamic)
- **HDF5 compression**: Reduces file size by ~50-70%

## Troubleshooting

### No snapshots group in HDF5 file

**Error**: `No snapshots found in HDF5 file`

**Solution**: Make sure `[output.snapshots]` is enabled in your config file:
```toml
[output.snapshots]
enabled = true
```

### Too many snapshots (file too large)

**Solution**: Increase snapshot intervals in config:
```toml
quasistatic_interval_years = 5.0   # Increase from 2.0
dynamic_interval = 0.2              # Increase from 0.1
```

### Not enough snapshots during earthquakes

**Solution**: Decrease dynamic interval or adjust velocity threshold:
```toml
dynamic_interval = 0.05             # More frequent sampling
velocity_threshold = 5.0e-4         # Lower threshold (0.5 mm/s)
```
