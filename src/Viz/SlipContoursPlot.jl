"""
    SlipContoursPlot

Plotting functions for cumulative slip contours (depth vs time).

These functions create space-time plots showing the evolution of slip
on the fault over time, similar to traditional earthquake cycle visualizations.
"""

using CairoMakie
using Printf


"""
    plot_slip_contours(h5_file::String;
                      output_file="slip_contours.png",
                      time_units="years",
                      figsize=(1000, 600),
                      max_depth=nothing,
                      colormap=:viridis,
                      levels=20)

Create contour plot of cumulative slip vs depth and time from snapshot data.

# Arguments
- `h5_file::String`: Path to HDF5 output file with snapshots
- `output_file::String`: Output filename (default: "slip_contours.png")
- `time_units::String`: "years" or "seconds" (default: "years")
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (1000, 600))
- `max_depth::Union{Float64,Nothing}`: Maximum depth to plot in km (default: nothing = all)
- `colormap::Symbol`: Colormap for contours (default: :viridis)
- `levels::Int`: Number of contour levels (default: 20)

# Example
```julia
using SEAS_SEME

# Plot slip contours
plot_slip_contours("data/strike_slip_2d/outputs/strike_slip_2d.h5")

# Custom settings
plot_slip_contours("data/strike_slip_2d/outputs/strike_slip_2d.h5",
                  time_units="years",
                  max_depth=30.0,
                  colormap=:thermal,
                  levels=30)
```
"""
function plot_slip_contours(h5_file::String;
                           output_file="slip_contours.png",
                           time_units="years",
                           figsize=(1000, 600),
                           max_depth=nothing,
                           colormap=:viridis,
                           levels=20)
    # Read snapshot data
    times, slip, slip_rate, stress, state, depths_km = read_snapshots(h5_file)

    # Convert time units
    yr2sec = 365.25 * 24 * 60 * 60
    if time_units == "years"
        times_plot = times ./ yr2sec
        time_label = "Time (years)"
    elseif time_units == "seconds"
        times_plot = times
        time_label = "Time (s)"
    else
        error("time_units must be 'years' or 'seconds'")
    end

    # Filter by depth if requested
    if !isnothing(max_depth)
        mask = depths_km .<= max_depth
        depths_km = depths_km[mask]
        slip = slip[mask, :]
    end

    # Create figure
    fig = Figure(resolution=figsize)

    ax = Axis(fig[1, 1],
             xlabel=time_label,
             ylabel="Depth (km)",
             yreversed=true,
             title="Cumulative Slip Evolution")

    # Create contour plot
    cf = contourf!(ax, times_plot, depths_km, slip,
                  levels=levels,
                  colormap=colormap)

    # Add colorbar
    Colorbar(fig[1, 2], cf, label="Cumulative Slip (m)")

    # Save figure
    save(output_file, fig)

    @info "Slip contours plot saved" file=output_file n_snapshots=length(times)

    return fig
end


"""
    plot_slip_rate_contours(h5_file::String;
                           output_file="slip_rate_contours.png",
                           time_units="years",
                           figsize=(1000, 600),
                           max_depth=nothing,
                           colormap=:hot,
                           levels=20,
                           log_scale=true)

Create contour plot of slip rate vs depth and time from snapshot data.

# Arguments
- `h5_file::String`: Path to HDF5 output file with snapshots
- `output_file::String`: Output filename (default: "slip_rate_contours.png")
- `time_units::String`: "years" or "seconds" (default: "years")
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (1000, 600))
- `max_depth::Union{Float64,Nothing}`: Maximum depth to plot in km (default: nothing = all)
- `colormap::Symbol`: Colormap for contours (default: :hot)
- `levels::Int`: Number of contour levels (default: 20)
- `log_scale::Bool`: Use logarithmic scale for slip rate (default: true)

# Example
```julia
plot_slip_rate_contours("data/strike_slip_2d/outputs/strike_slip_2d.h5",
                       log_scale=true,
                       colormap=:hot)
```
"""
function plot_slip_rate_contours(h5_file::String;
                                output_file="slip_rate_contours.png",
                                time_units="years",
                                figsize=(1000, 600),
                                max_depth=nothing,
                                colormap=:hot,
                                levels=20,
                                log_scale=true)
    # Read snapshot data
    times, slip, slip_rate, stress, state, depths_km = read_snapshots(h5_file)

    # Convert time units
    yr2sec = 365.25 * 24 * 60 * 60
    if time_units == "years"
        times_plot = times ./ yr2sec
        time_label = "Time (years)"
    elseif time_units == "seconds"
        times_plot = times
        time_label = "Time (s)"
    else
        error("time_units must be 'years' or 'seconds'")
    end

    # Filter by depth if requested
    if !isnothing(max_depth)
        mask = depths_km .<= max_depth
        depths_km = depths_km[mask]
        slip_rate = slip_rate[mask, :]
    end

    # Apply log scale if requested
    if log_scale
        # Avoid log(0) by setting minimum value
        slip_rate_plot = log10.(max.(slip_rate, 1e-12))
        colorbar_label = "log₁₀(Slip Rate) [m/s]"
    else
        slip_rate_plot = slip_rate
        colorbar_label = "Slip Rate (m/s)"
    end

    # Create figure
    fig = Figure(resolution=figsize)

    ax = Axis(fig[1, 1],
             xlabel=time_label,
             ylabel="Depth (km)",
             yreversed=true,
             title="Slip Rate Evolution")

    # Create contour plot
    cf = contourf!(ax, times_plot, depths_km, slip_rate_plot,
                  levels=levels,
                  colormap=colormap)

    # Add colorbar
    Colorbar(fig[1, 2], cf, label=colorbar_label)

    # Save figure
    save(output_file, fig)

    @info "Slip rate contours plot saved" file=output_file n_snapshots=length(times)

    return fig
end


"""
    plot_stress_contours(h5_file::String;
                        output_file="stress_contours.png",
                        time_units="years",
                        figsize=(1000, 600),
                        max_depth=nothing,
                        colormap=:balance,
                        levels=20)

Create contour plot of shear stress vs depth and time from snapshot data.

# Arguments
- `h5_file::String`: Path to HDF5 output file with snapshots
- `output_file::String`: Output filename (default: "stress_contours.png")
- `time_units::String`: "years" or "seconds" (default: "years")
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (1000, 600))
- `max_depth::Union{Float64,Nothing}`: Maximum depth to plot in km (default: nothing = all)
- `colormap::Symbol`: Colormap for contours (default: :balance)
- `levels::Int`: Number of contour levels (default: 20)

# Example
```julia
plot_stress_contours("data/strike_slip_2d/outputs/strike_slip_2d.h5")
```
"""
function plot_stress_contours(h5_file::String;
                             output_file="stress_contours.png",
                             time_units="years",
                             figsize=(1000, 600),
                             max_depth=nothing,
                             colormap=:balance,
                             levels=20)
    # Read snapshot data
    times, slip, slip_rate, stress, state, depths_km = read_snapshots(h5_file)

    # Convert time units
    yr2sec = 365.25 * 24 * 60 * 60
    if time_units == "years"
        times_plot = times ./ yr2sec
        time_label = "Time (years)"
    elseif time_units == "seconds"
        times_plot = times
        time_label = "Time (s)"
    else
        error("time_units must be 'years' or 'seconds'")
    end

    # Filter by depth if requested
    if !isnothing(max_depth)
        mask = depths_km .<= max_depth
        depths_km = depths_km[mask]
        stress = stress[mask, :]
    end

    # Create figure
    fig = Figure(resolution=figsize)

    ax = Axis(fig[1, 1],
             xlabel=time_label,
             ylabel="Depth (km)",
             yreversed=true,
             title="Shear Stress Evolution")

    # Create contour plot
    cf = contourf!(ax, times_plot, depths_km, stress,
                  levels=levels,
                  colormap=colormap)

    # Add colorbar
    Colorbar(fig[1, 2], cf, label="Shear Stress (MPa)")

    # Save figure
    save(output_file, fig)

    @info "Stress contours plot saved" file=output_file n_snapshots=length(times)

    return fig
end


"""
    plot_all_contours(h5_file::String;
                     output_dir=".",
                     time_units="years",
                     max_depth=nothing)

Generate all contour plots (slip, slip rate, stress) in one call.

# Arguments
- `h5_file::String`: Path to HDF5 output file with snapshots
- `output_dir::String`: Directory for output files (default: current directory)
- `time_units::String`: "years" or "seconds" (default: "years")
- `max_depth::Union{Float64,Nothing}`: Maximum depth to plot in km (default: nothing = all)

# Example
```julia
# Generate all plots
plot_all_contours("data/strike_slip_2d/outputs/strike_slip_2d.h5",
                 output_dir="data/strike_slip_2d/plots",
                 max_depth=30.0)
```
"""
function plot_all_contours(h5_file::String;
                          output_dir=".",
                          time_units="years",
                          max_depth=nothing)
    mkpath(output_dir)

    @info "Generating contour plots..." output_dir

    # Slip contours
    plot_slip_contours(h5_file,
                      output_file=joinpath(output_dir, "slip_contours.png"),
                      time_units=time_units,
                      max_depth=max_depth)

    # Slip rate contours
    plot_slip_rate_contours(h5_file,
                           output_file=joinpath(output_dir, "slip_rate_contours.png"),
                           time_units=time_units,
                           max_depth=max_depth)

    # Stress contours
    plot_stress_contours(h5_file,
                        output_file=joinpath(output_dir, "stress_contours.png"),
                        time_units=time_units,
                        max_depth=max_depth)

    @info "All contour plots generated" directory=output_dir

    return nothing
end
