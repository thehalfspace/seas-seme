"""
    Viz

Visualization module for SEAS-SEME simulation results.

This module provides modular, reusable plotting functions for:
- Parameter distributions (friction, stress, initial conditions)
- Maximum slip rate time series
- Earthquake cycle heatmaps
- Cumulative slip profiles

# Usage

```julia
using SEAS_SEME.Viz

# Plot initial conditions
plot_initial_conditions("data/strike_slip_benchmark/params")

# Plot maximum slip rate
plot_vfmax("data/strike_slip_benchmark/outputs/strike_slip_benchmark.h5")

# Plot earthquake cycle heatmap
plot_earthquake_cycle("data/strike_slip_benchmark/outputs/strike_slip_benchmark.h5",
                     depth_range=(0.0, 20.0))

# Plot cumulative slip
plot_cumulative_slip("data/strike_slip_benchmark/outputs/strike_slip_benchmark.h5")
```

# Exported Functions

## Data Readers
- `read_timeseries_data`: Read global time series (time, Vfmax)
- `read_fault_geometry`: Read fault geometry
- `read_depth_timeseries`: Read time series at specific depth
- `read_initial_conditions`: Read friction parameters and stresses
- `read_all_fault_timeseries`: Read slip rate matrix for heatmaps
- `get_available_depths`: Get list of available depths
- `read_metadata`: Read simulation metadata
- `read_snapshots`: Read all snapshot data (slip, slip rate, stress, state)
- `read_snapshot_at_index`: Read single snapshot at specific index
- `get_snapshot_config`: Get snapshot configuration from file

## Parameter Plots
- `plot_initial_conditions`: Plot friction + stress with dual axes
- `plot_friction_parameters`: Plot a, b, Lc in separate panels
- `plot_stress_distribution`: Plot normal and shear stress
- `plot_ab_difference`: Plot (a-b) showing VW/VS regions

## Vfmax Plots
- `plot_vfmax`: Plot maximum slip rate vs time
- `plot_vfmax_comparison`: Compare Vfmax from multiple simulations
- `plot_vfmax_with_events`: Plot Vfmax with earthquake events highlighted

## Earthquake Cycle Plots
- `plot_earthquake_cycle`: Slip rate heatmap (depth vs time)
- `plot_earthquake_cycle_with_vfmax`: Combined heatmap + Vfmax plot
- `plot_depth_slice`: Slip rate time series at specific depth

## Cumulative Slip Plots
- `plot_cumulative_slip`: Slip profiles at multiple times
- `plot_slip_at_times`: Slip profiles at specified times
- `plot_slip_evolution`: Cumulative slip vs time at specific depth
- `plot_slip_deficit`: Slip deficit (plate motion - slip) vs time

## Slip Contour Plots (from Snapshots)
- `plot_slip_contours`: Contour plot of cumulative slip (depth vs time)
- `plot_slip_rate_contours`: Contour plot of slip rate (depth vs time)
- `plot_stress_contours`: Contour plot of shear stress (depth vs time)
- `plot_all_contours`: Generate all contour plots in one call
"""
module Viz

using CairoMakie
using HDF5
using DelimitedFiles
using Printf

# Include all plotting modules
include("DataReader.jl")
include("ParametersPlot.jl")
include("VfmaxPlot.jl")
include("EarthquakeCyclePlot.jl")
include("CumulativeSlipPlot.jl")
include("SlipContoursPlot.jl")

# Export data reading functions
export read_timeseries_data
export read_fault_geometry
export read_depth_timeseries
export read_initial_conditions
export read_all_fault_timeseries
export get_available_depths
export read_metadata
export read_snapshots
export read_snapshot_at_index
export get_snapshot_config

# Export parameter plotting functions
export plot_initial_conditions
export plot_friction_parameters
export plot_stress_distribution
export plot_ab_difference

# Export Vfmax plotting functions
export plot_vfmax
export plot_vfmax_comparison
export plot_vfmax_with_events

# Export earthquake cycle plotting functions
export plot_earthquake_cycle
export plot_earthquake_cycle_with_vfmax
export plot_depth_slice

# Export cumulative slip plotting functions
export plot_cumulative_slip
export plot_slip_at_times
export plot_slip_evolution
export plot_slip_deficit

# Export slip contour plotting functions
export plot_slip_contours
export plot_slip_rate_contours
export plot_stress_contours
export plot_all_contours

end  # module Viz
