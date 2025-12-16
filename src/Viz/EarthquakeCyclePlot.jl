"""
    EarthquakeCyclePlot

Plotting functions for earthquake cycle visualization (slip rate heatmap).
"""

using CairoMakie
using Printf


"""
    plot_earthquake_cycle(h5_file::String;
                         output_file="eq_cycle.png",
                         depth_range=nothing,
                         time_units="timesteps",
                         figsize=(800, 600),
                         log_scale=true,
                         vmin=1e-9,
                         vmax=1e0,
                         colormap=:inferno,
                         interpolate=true)

Create slip rate heatmap (depth vs time) to visualize earthquake cycles.

# Arguments
- `h5_file::String`: Path to HDF5 output file
- `output_file::String`: Output filename (default: "eq_cycle.png")
- `depth_range::Union{Tuple{Float64,Float64},Nothing}`: Depth range in km (default: nothing = all)
- `time_units::String`: "timesteps" or "years" (default: "timesteps")
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (800, 600))
- `log_scale::Bool`: Use log scale for color (default: true)
- `vmin::Float64`: Minimum slip rate for colorbar in m/s (default: 1e-9)
- `vmax::Float64`: Maximum slip rate for colorbar in m/s (default: 1e0)
- `colormap::Symbol`: Colormap name (default: :inferno)
- `interpolate::Bool`: Interpolate heatmap (default: true)

# Example
```julia
plot_earthquake_cycle("outputs/strike_slip_benchmark.h5")
plot_earthquake_cycle("outputs/strike_slip_benchmark.h5", depth_range=(0.0, 20.0))
plot_earthquake_cycle("outputs/strike_slip_benchmark.h5", time_units="years")
```
"""
function plot_earthquake_cycle(h5_file::String;
                              output_file="eq_cycle.png",
                              depth_range=nothing,
                              time_units="timesteps",
                              figsize=(800, 600),
                              log_scale=true,
                              vmin=1e-9,
                              vmax=1e0,
                              colormap=:inferno,
                              interpolate=true)
    # Read data
    times, slip_rate_matrix, depths_km = read_all_fault_timeseries(h5_file)

    # Filter by depth range if specified
    if !isnothing(depth_range)
        mask = (depths_km .>= depth_range[1]) .& (depths_km .<= depth_range[2])
        slip_rate_matrix = slip_rate_matrix[mask, :]
        depths_km = depths_km[mask]
    end

    # Prepare time axis
    yr2sec = 365.25 * 24 * 60 * 60
    if time_units == "years"
        time_values = times ./ yr2sec
        time_label = "Time (years)"
    elseif time_units == "timesteps"
        time_values = 1:length(times)
        time_label = "Timesteps"
    else
        error("time_units must be 'years' or 'timesteps'")
    end

    # Create figure
    fig = Figure(resolution=figsize)

    ax = Axis(fig[1, 1],
             xlabel=time_label,
             ylabel="Depth (km)",
             yreversed=true,
             title="Earthquake Cycle (Slip Rate)")

    # Create heatmap
    if log_scale
        # Use log normalization
        hm = heatmap!(ax, time_values, depths_km, slip_rate_matrix,
                     colormap=colormap,
                     colorrange=(vmin, vmax),
                     colorscale=log10,
                     interpolate=interpolate)
    else
        hm = heatmap!(ax, time_values, depths_km, slip_rate_matrix,
                     colormap=colormap,
                     colorrange=(vmin, vmax),
                     interpolate=interpolate)
    end

    # Add colorbar
    Colorbar(fig[1, 2], hm, label="Slip Rate (m/s)")

    # Save figure
    save(output_file, fig)

    @info "Earthquake cycle plot saved" file=output_file

    return fig
end


"""
    plot_earthquake_cycle_with_vfmax(h5_file::String;
                                    output_file="eq_cycle_with_vfmax.png",
                                    depth_range=nothing,
                                    figsize=(800, 800),
                                    log_scale=true,
                                    vmin=1e-9,
                                    vmax=1e0,
                                    colormap=:inferno)

Create combined plot: slip rate heatmap + Vfmax time series.

Creates a two-panel vertical plot with:
- Top panel: Vfmax vs time
- Bottom panel: Slip rate heatmap (depth vs time)

# Arguments
- `h5_file::String`: Path to HDF5 output file
- `output_file::String`: Output filename (default: "eq_cycle_with_vfmax.png")
- `depth_range::Union{Tuple{Float64,Float64},Nothing}`: Depth range in km (default: nothing = all)
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (800, 800))
- `log_scale::Bool`: Use log scale for color (default: true)
- `vmin::Float64`: Minimum slip rate for colorbar in m/s (default: 1e-9)
- `vmax::Float64`: Maximum slip rate for colorbar in m/s (default: 1e0)
- `colormap::Symbol`: Colormap name (default: :inferno)

# Example
```julia
plot_earthquake_cycle_with_vfmax("outputs/strike_slip_benchmark.h5")
```
"""
function plot_earthquake_cycle_with_vfmax(h5_file::String;
                                         output_file="eq_cycle_with_vfmax.png",
                                         depth_range=nothing,
                                         figsize=(800, 800),
                                         log_scale=true,
                                         vmin=1e-9,
                                         vmax=1e0,
                                         colormap=:inferno)
    # Read data
    times_global, Vfmax = read_timeseries_data(h5_file)
    times, slip_rate_matrix, depths_km = read_all_fault_timeseries(h5_file)

    # Convert to years
    yr2sec = 365.25 * 24 * 60 * 60
    times_years = times ./ yr2sec
    times_global_years = times_global ./ yr2sec

    # Filter by depth range if specified
    if !isnothing(depth_range)
        mask = (depths_km .>= depth_range[1]) .& (depths_km .<= depth_range[2])
        slip_rate_matrix = slip_rate_matrix[mask, :]
        depths_km = depths_km[mask]
    end

    # Create figure with 2 rows
    fig = Figure(resolution=figsize)

    # Top panel: Vfmax
    ax1 = Axis(fig[1, 1],
              ylabel="Vfmax (m/s)",
              yscale=log10,
              title="Maximum Fault Slip Rate")

    lines!(ax1, times_global_years, Vfmax, color=:blue, linewidth=2)

    # Hide x-axis labels for top panel
    hidexdecorations!(ax1, grid=false)

    # Bottom panel: Heatmap
    ax2 = Axis(fig[2, 1],
              xlabel="Time (years)",
              ylabel="Depth (km)",
              yreversed=true)

    if log_scale
        hm = heatmap!(ax2, times_years, depths_km, slip_rate_matrix,
                     colormap=colormap,
                     colorrange=(vmin, vmax),
                     colorscale=log10,
                     interpolate=true)
    else
        hm = heatmap!(ax2, times_years, depths_km, slip_rate_matrix,
                     colormap=colormap,
                     colorrange=(vmin, vmax),
                     interpolate=true)
    end

    # Add colorbar
    Colorbar(fig[2, 2], hm, label="Slip Rate (m/s)")

    # Link x-axes
    linkxaxes!(ax1, ax2)

    # Save figure
    save(output_file, fig)

    @info "Earthquake cycle with Vfmax plot saved" file=output_file

    return fig
end


"""
    plot_depth_slice(h5_file::String, depth::Float64;
                    output_file="depth_slice.png",
                    time_units="years",
                    figsize=(800, 400),
                    log_scale=true)

Plot slip rate time series at a specific depth.

# Arguments
- `h5_file::String`: Path to HDF5 output file
- `depth::Float64`: Depth in km
- `output_file::String`: Output filename (default: "depth_slice.png")
- `time_units::String`: "years" or "seconds" (default: "years")
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (800, 400))
- `log_scale::Bool`: Use log scale for y-axis (default: true)

# Example
```julia
plot_depth_slice("outputs/strike_slip_benchmark.h5", 10.0)
```
"""
function plot_depth_slice(h5_file::String, depth::Float64;
                         output_file="depth_slice.png",
                         time_units="years",
                         figsize=(800, 400),
                         log_scale=true)
    # Read data for specific depth
    times, slip, slip_rate, stress, state = read_depth_timeseries(h5_file, depth)

    # Convert time units
    yr2sec = 365.25 * 24 * 60 * 60
    if time_units == "years"
        times = times ./ yr2sec
        time_label = "Time (years)"
    elseif time_units == "seconds"
        time_label = "Time (s)"
    else
        error("time_units must be 'years' or 'seconds'")
    end

    # Create figure
    fig = Figure(resolution=figsize)

    ax = Axis(fig[1, 1],
             xlabel=time_label,
             ylabel="Slip Rate (m/s)",
             title="Slip Rate at $(depth) km depth")

    if log_scale
        ax.yscale = log10
    end

    lines!(ax, times, slip_rate, color=:blue, linewidth=2)

    # Save figure
    save(output_file, fig)

    @info "Depth slice plot saved" file=output_file depth=depth

    return fig
end
