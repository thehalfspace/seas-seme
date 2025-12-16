"""
Example script demonstrating the Viz module usage.

This script shows how to use the SEAS_SEME.Viz module to create
various plots from simulation results.
"""

using SEAS_SEME.Viz

# Set paths
simulation_dir = "data/strike_slip_benchmark"
params_dir = joinpath(simulation_dir, "params")
h5_file = joinpath(simulation_dir, "outputs/strike_slip_benchmark.h5")
output_dir = "plots/strike_slip_benchmark"

# Create output directory
mkpath(output_dir)

println("Generating visualization plots...")

# ========================================
# 1. Parameter Plots
# ========================================

println("\n1. Plotting initial conditions and parameters...")

# Plot initial conditions (friction + stress with dual axes)
plot_initial_conditions(params_dir,
                       output_file=joinpath(output_dir, "initial_conditions.png"))

# Plot friction parameters in separate panels
plot_friction_parameters(params_dir,
                        output_file=joinpath(output_dir, "friction_parameters.png"))

# Plot stress distribution
plot_stress_distribution(params_dir,
                        output_file=joinpath(output_dir, "stress_distribution.png"))

# Plot (a-b) to show VW/VS regions
plot_ab_difference(params_dir,
                  output_file=joinpath(output_dir, "ab_difference.png"))

# ========================================
# 2. Vfmax Plots
# ========================================

println("\n2. Plotting maximum slip rate (Vfmax)...")

# Basic Vfmax plot
plot_vfmax(h5_file,
          output_file=joinpath(output_dir, "vfmax.png"),
          time_units="years",
          log_scale=true)

# Vfmax with earthquake events highlighted
plot_vfmax_with_events(h5_file,
                      output_file=joinpath(output_dir, "vfmax_events.png"),
                      event_threshold=1e-3)

# ========================================
# 3. Earthquake Cycle Plots
# ========================================

println("\n3. Plotting earthquake cycle heatmaps...")

# Basic earthquake cycle heatmap
plot_earthquake_cycle(h5_file,
                     output_file=joinpath(output_dir, "eq_cycle.png"),
                     depth_range=(0.0, 20.0),
                     time_units="timesteps",
                     log_scale=true)

# Combined plot: Vfmax + heatmap
plot_earthquake_cycle_with_vfmax(h5_file,
                                output_file=joinpath(output_dir, "eq_cycle_with_vfmax.png"),
                                depth_range=(0.0, 20.0))

# Plot slip rate at a specific depth
plot_depth_slice(h5_file, 10.0,
                output_file=joinpath(output_dir, "depth_10km_sliprate.png"),
                time_units="years")

# ========================================
# 4. Cumulative Slip Plots
# ========================================

println("\n4. Plotting cumulative slip...")

# Cumulative slip profiles at multiple times
plot_cumulative_slip(h5_file,
                    output_file=joinpath(output_dir, "cumulative_slip.png"),
                    dynamic_interval=0.1,
                    quasistatic_interval=2.0)

# Slip profiles at specific times
plot_slip_at_times(h5_file, [10.0, 50.0, 100.0, 200.0],
                  output_file=joinpath(output_dir, "slip_at_times.png"),
                  time_units="years")

# Cumulative slip evolution at a specific depth
plot_slip_evolution(h5_file, 10.0,
                   output_file=joinpath(output_dir, "slip_evolution_10km.png"),
                   time_units="years")

# Slip deficit (requires plate velocity)
# Assuming Vpl = 1e-9 m/s (approximately 3.15 cm/yr)
plot_slip_deficit(h5_file, 1e-9,
                 output_file=joinpath(output_dir, "slip_deficit_10km.png"),
                 depth=10.0)

println("\n✓ All plots generated successfully!")
println("Output directory: $output_dir")

# ========================================
# 5. Advanced: Custom Plots
# ========================================

println("\n5. Demonstrating data access for custom plots...")

# You can also read data directly and create custom plots
times, Vfmax = read_timeseries_data(h5_file)
println("  - Read $(length(times)) timesteps")

depths = get_available_depths(h5_file)
println("  - Available depths: $(length(depths)) locations")

meta = read_metadata(h5_file)
println("  - Simulation: $(meta["simulation_name"])")
println("  - Total time: $(meta["total_time"]) s")

println("\nExample complete!")
