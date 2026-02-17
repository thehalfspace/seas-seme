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
                      figsize=(800, 600),
                      max_depth=nothing,
                      velocity_threshold=1e-3,
                      dynamic_step=1,
                      quasistatic_step=1)

Plot cumulative slip profiles from snapshot data with seismic/interseismic coloring.

- Y-axis: Depth (km), surface at top
- X-axis: Cumulative slip (m)
- Red lines: seismic (dynamic) snapshots (max slip rate >= threshold)
- Blue lines: interseismic (quasistatic) snapshots

# Arguments
- `h5_file::String`: Path to HDF5 output file with snapshots
- `output_file::String`: Output filename (default: "slip_contours.png")
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (800, 600))
- `max_depth::Union{Float64,Nothing}`: Maximum depth to plot in km (default: nothing = all)
- `velocity_threshold::Float64`: Slip rate threshold for seismic classification (m/s, default: 1e-3)
- `dynamic_step::Int`: Plot every Nth dynamic snapshot (default: 1 = all)
- `quasistatic_step::Int`: Plot every Nth quasistatic snapshot (default: 1 = all)

# Example
```julia
plot_slip_contours("data/plane_strain_2d/outputs/plane_strain_2d.h5",
                  max_depth=20.0, dynamic_step=10)
```
"""
function plot_slip_contours(h5_file::String;
                           output_file="slip_contours.png",
                           figsize=(800, 600),
                           max_depth=nothing,
                           velocity_threshold=1e-3,
                           dynamic_step=1,
                           quasistatic_step=1)
    # Read snapshot data: slip is [n_fault x n_snapshots]
    times, slip, slip_rate, stress, state, depths_km = read_snapshots(h5_file)

    # Filter by depth if requested
    if !isnothing(max_depth)
        mask = depths_km .<= max_depth
        depths_km = depths_km[mask]
        slip = slip[mask, :]
        slip_rate = slip_rate[mask, :]
    end

    # Shift depths so free surface = 0 (surface is max value in coordinate system)
    #surface_depth = maximum(depths_km)
    #depths_km = depths_km .- surface_depth  # e.g., -40 → 0, -20 → 20

    # Classify each snapshot by max slip rate on fault
    max_sr = vec(maximum(slip_rate, dims=1))
    is_dynamic = max_sr .>= velocity_threshold

    dynamic_indices = findall(is_dynamic)
    quasistatic_indices = findall(.!is_dynamic)

    # Subsample
    dynamic_indices = dynamic_indices[1:dynamic_step:end]
    quasistatic_indices = quasistatic_indices[1:quasistatic_step:end]

    # Create figure
    fig = Figure(size=figsize)

    ax = Axis(fig[1, 1],
             xlabel="Cumulative Slip (m)",
             ylabel="Depth (km)",
             yreversed=true,
             title="Cumulative Slip Profiles")

    # Plot interseismic first (behind)
    for idx in quasistatic_indices
        lines!(ax, slip[:, idx], depths_km,
              color=(:royalblue, 0.6), linewidth=0.8)
    end

    # Plot seismic on top
    for idx in dynamic_indices
        lines!(ax, slip[:, idx], depths_km,
              color=(:chocolate, 0.5), linewidth=1.2)
    end

    # Legend
    lines!(ax, [NaN], [NaN], color=:chocolate, linewidth=1.5, label="Seismic")
    lines!(ax, [NaN], [NaN], color=:royalblue, linewidth=1.0, label="Interseismic")
    axislegend(ax, position=:rb)

    # Save figure
    save(output_file, fig)

    @info "Slip contours plot saved" file=output_file n_dynamic=length(dynamic_indices) n_quasistatic=length(quasistatic_indices)

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
