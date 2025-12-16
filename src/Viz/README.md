# Viz Module

Modular visualization tools for SEAS-SEME simulation results.

## Overview

The Viz module provides publication-quality plotting functions for earthquake cycle simulation results using CairoMakie. All plots are organized into logical categories:

- **Parameter Plots**: Friction parameters, stress distributions, initial conditions
- **Vfmax Plots**: Maximum slip rate time series
- **Earthquake Cycle Plots**: Slip rate heatmaps (depth vs time)
- **Cumulative Slip Plots**: Slip profiles and evolution

## Quick Start

```julia
using SEAS_SEME.Viz

# Set paths
params_dir = "data/strike_slip_benchmark/params"
h5_file = "data/strike_slip_benchmark/outputs/strike_slip_benchmark.h5"

# Plot initial conditions
plot_initial_conditions(params_dir, output_file="initial_conditions.png")

# Plot maximum slip rate
plot_vfmax(h5_file, output_file="vfmax.png", time_units="years")

# Plot earthquake cycle heatmap
plot_earthquake_cycle(h5_file, output_file="eq_cycle.png", depth_range=(0.0, 20.0))

# Plot cumulative slip profiles
plot_cumulative_slip(h5_file, output_file="cumulative_slip.png")
```

## Module Structure

```
src/Viz/
├── DataReader.jl          # HDF5 and parameter file reading utilities
├── ParametersPlot.jl      # Friction, stress, and initial condition plots
├── VfmaxPlot.jl          # Maximum slip rate plots
├── EarthquakeCyclePlot.jl # Slip rate heatmaps
├── CumulativeSlipPlot.jl  # Cumulative slip profiles
└── Viz.jl                # Main module interface
```

## Data Readers

### `read_timeseries_data(h5_file)`
Read global time series (time, Vfmax) from HDF5 file.

### `read_fault_geometry(h5_file)`
Read fault coordinates and depths from HDF5 file.

### `read_depth_timeseries(h5_file, depth)`
Read time series data at a specific depth (slip, slip_rate, stress, state).

### `read_initial_conditions(params_dir)`
Read friction parameters (a, b, Lc) and stresses (σ_n, τ) from params directory.

### `read_all_fault_timeseries(h5_file)`
Read slip rate matrix for all fault nodes (for heatmaps).

### `get_available_depths(h5_file)`
Get list of available depth locations in HDF5 file.

### `read_metadata(h5_file)`
Read simulation metadata from HDF5 file.

## Parameter Plots

### `plot_initial_conditions(params_dir; kwargs...)`
Plot friction parameters and initial stresses with dual axes.

**Arguments:**
- `params_dir::String`: Path to params directory
- `output_file`: Output filename (default: "initial_conditions.png")
- `figsize`: Figure size in pixels (default: (800, 600))
- `max_depth`: Maximum depth to plot in km (default: nothing = all)

### `plot_friction_parameters(params_dir; kwargs...)`
Plot friction parameters (a, b, Lc) in separate panels.

### `plot_stress_distribution(params_dir; kwargs...)`
Plot normal and shear stress vs depth.

### `plot_ab_difference(params_dir; kwargs...)`
Plot (a-b) to show velocity-weakening/strengthening regions.

## Vfmax Plots

### `plot_vfmax(h5_file; kwargs...)`
Plot maximum fault slip rate vs time.

**Arguments:**
- `h5_file::String`: Path to HDF5 file
- `output_file`: Output filename (default: "vfmax.png")
- `time_units`: "years" or "seconds" (default: "years")
- `figsize`: Figure size in pixels (default: (800, 400))
- `log_scale`: Use log scale for y-axis (default: true)
- `time_range`: Time range to plot (default: nothing = all)

### `plot_vfmax_comparison(h5_files, labels; kwargs...)`
Compare Vfmax from multiple simulations.

### `plot_vfmax_with_events(h5_file; kwargs...)`
Plot Vfmax with earthquake events highlighted.

## Earthquake Cycle Plots

### `plot_earthquake_cycle(h5_file; kwargs...)`
Create slip rate heatmap (depth vs time).

**Arguments:**
- `h5_file::String`: Path to HDF5 file
- `output_file`: Output filename (default: "eq_cycle.png")
- `depth_range`: Depth range in km (default: nothing = all)
- `time_units`: "timesteps" or "years" (default: "timesteps")
- `figsize`: Figure size in pixels (default: (800, 600))
- `log_scale`: Use log scale for colorbar (default: true)
- `vmin`: Minimum slip rate for colorbar (default: 1e-9)
- `vmax`: Maximum slip rate for colorbar (default: 1e0)
- `colormap`: Colormap name (default: :inferno)

### `plot_earthquake_cycle_with_vfmax(h5_file; kwargs...)`
Combined plot: Vfmax + slip rate heatmap.

### `plot_depth_slice(h5_file, depth; kwargs...)`
Plot slip rate time series at a specific depth.

## Cumulative Slip Plots

### `plot_cumulative_slip(h5_file; kwargs...)`
Plot cumulative slip profiles at multiple times.

**Arguments:**
- `h5_file::String`: Path to HDF5 file
- `output_file`: Output filename (default: "cumulative_slip.png")
- `dynamic_interval`: Sampling interval during events in seconds (default: 0.1)
- `quasistatic_interval`: Sampling interval during quasistatic in years (default: 2.0)
- `figsize`: Figure size in pixels (default: (800, 600))
- `velocity_threshold`: Threshold for identifying events (default: 1e-3 m/s)

### `plot_slip_at_times(h5_file, times; kwargs...)`
Plot slip profiles at specified times.

### `plot_slip_evolution(h5_file, depth; kwargs...)`
Plot cumulative slip vs time at a specific depth.

### `plot_slip_deficit(h5_file, plate_velocity; kwargs...)`
Plot slip deficit (plate motion - slip) vs time.

## Example Usage

See `examples/viz_example.jl` for a complete example demonstrating all plotting functions.

```bash
julia examples/viz_example.jl
```

## Design Principles

1. **Modular**: Each plot type in its own file
2. **Separation of concerns**: Data reading separate from plotting
3. **Reusable**: Functions accept file paths, return figures
4. **Configurable**: Key parameters (figsize, colors, intervals) as optional args
5. **Self-contained**: Each function can be used independently
6. **Modern stack**: CairoMakie for plotting, HDF5.jl for data

## Notes

- All plotting uses CairoMakie (modern, fast, publication-quality)
- Data reading uses HDF5.jl (efficient binary format)
- Follows SEAS-SEME modular architecture pattern
- Outputs saved as PNG (can be extended to PDF, SVG)
- All depth values are in kilometers (positive down)
- All stress values are in MPa
- All time values can be in years or seconds (configurable)
