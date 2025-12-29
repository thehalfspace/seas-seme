"""
    CumulativeSlipPlot

Plotting functions for cumulative slip profiles vs depth.
"""

using CairoMakie
using Printf


"""
    plot_cumulative_slip(h5_file::String;
                        output_file="cumulative_slip.png",
                        dynamic_interval=0.1,
                        quasistatic_interval=2.0,
                        figsize=(800, 600),
                        max_depth=nothing,
                        velocity_threshold=1e-3)

Plot cumulative slip profiles at different times.

Samples slip profiles:
- Every `dynamic_interval` seconds during dynamic events (Vfmax > threshold)
- Every `quasistatic_interval` years during quasistatic periods

# Arguments
- `h5_file::String`: Path to HDF5 output file
- `output_file::String`: Output filename (default: "cumulative_slip.png")
- `dynamic_interval::Float64`: Sampling interval during events in seconds (default: 0.1)
- `quasistatic_interval::Float64`: Sampling interval during quasistatic in years (default: 2.0)
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (800, 600))
- `max_depth::Union{Float64,Nothing}`: Maximum depth to plot in km (default: nothing = all)
- `velocity_threshold::Float64`: Threshold for identifying dynamic events in m/s (default: 1e-3)

# Example
```julia
plot_cumulative_slip("outputs/strike_slip_benchmark.h5")
plot_cumulative_slip("outputs/strike_slip_benchmark.h5", dynamic_interval=0.05, quasistatic_interval=1.0)
```
"""
function plot_cumulative_slip(h5_file::String;
                             output_file="cumulative_slip.png",
                             dynamic_interval=0.1,
                             quasistatic_interval=2.0,
                             figsize=(800, 600),
                             max_depth=nothing,
                             velocity_threshold=1e-3)
    # Read global time series to identify dynamic periods
    times_global, Vfmax, solver_mode = read_timeseries_data(h5_file)

    # Get available depths
    available_depths = get_available_depths(h5_file)

    # Read first depth to get all time steps
    _, depths_km = read_fault_geometry(h5_file)

    # Filter by max depth
    if !isnothing(max_depth)
        mask = depths_km .<= max_depth
        depths_km = depths_km[mask]
    end

    yr2sec = 365.25 * 24 * 60 * 60

    # Use solver_mode if available, otherwise fall back to velocity threshold
    use_solver_mode = !isempty(solver_mode)

    # Identify sampling indices
    sampling_indices = Int[]
    last_dynamic_time = -Inf
    last_quasistatic_time = -Inf

    for (i, t) in enumerate(times_global)
        # Determine if dynamic based on solver_mode or velocity
        is_dynamic = if use_solver_mode
            solver_mode[i] == 1
        else
            Vfmax[i] >= velocity_threshold
        end

        if is_dynamic
            # Dynamic event
            if t - last_dynamic_time >= dynamic_interval
                push!(sampling_indices, i)
                last_dynamic_time = t
            end
        else
            # Quasistatic period
            if t - last_quasistatic_time >= quasistatic_interval * yr2sec
                push!(sampling_indices, i)
                last_quasistatic_time = t
            end
        end
    end

    # Create figure
    fig = Figure(resolution=figsize)

    ax = Axis(fig[1, 1],
             xlabel="Accumulated Slip (m)",
             ylabel="Depth (km)",
             yreversed=true,
             title="Cumulative Slip Profiles")

    # Plot slip profiles for sampled time steps
    for (count, idx) in enumerate(sampling_indices)
        # Read slip at all depths for this time step
        slip_profile = zeros(length(available_depths))

        for (j, depth) in enumerate(available_depths)
            time, slip, _, _, _ = read_depth_timeseries(h5_file, depth)
            slip_profile[j] = slip[idx]
        end

        # Determine color based on solver_mode or velocity
        is_dynamic = if use_solver_mode
            solver_mode[idx] == 1
        else
            Vfmax[idx] >= velocity_threshold
        end

        if is_dynamic
            color = :chocolate  # Dynamic
            alpha = 0.5
            linewidth = 1.5
        else
            color = :royalblue  # Quasistatic
            alpha = 0.7
            linewidth = 1.0
        end

        lines!(ax, slip_profile, available_depths,
              color=(color, alpha), linewidth=linewidth)
    end

    # Add legend
    lines!(ax, [NaN], [NaN], color=:chocolate, linewidth=1.5, label="Dynamic")
    lines!(ax, [NaN], [NaN], color=:royalblue, linewidth=1.0, label="Quasistatic")
    axislegend(ax, position=:rb)

    # Save figure
    save(output_file, fig)

    @info "Cumulative slip plot saved" file=output_file n_profiles=length(sampling_indices)

    return fig
end


"""
    plot_slip_at_times(h5_file::String, times::Vector{Float64};
                      output_file="slip_at_times.png",
                      time_units="years",
                      figsize=(800, 600),
                      max_depth=nothing)

Plot cumulative slip profiles at specified times.

# Arguments
- `h5_file::String`: Path to HDF5 output file
- `times::Vector{Float64}`: Times to plot (in units specified by time_units)
- `output_file::String`: Output filename (default: "slip_at_times.png")
- `time_units::String`: "years" or "seconds" (default: "years")
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (800, 600))
- `max_depth::Union{Float64,Nothing}`: Maximum depth to plot in km (default: nothing = all)

# Example
```julia
plot_slip_at_times("outputs/strike_slip_benchmark.h5", [10.0, 50.0, 100.0, 200.0])
```
"""
function plot_slip_at_times(h5_file::String, times::Vector{Float64};
                           output_file="slip_at_times.png",
                           time_units="years",
                           figsize=(800, 600),
                           max_depth=nothing)
    # Get available depths
    available_depths = get_available_depths(h5_file)

    # Read time vector
    times_global, _, _ = read_timeseries_data(h5_file)

    # Convert requested times to seconds if needed
    yr2sec = 365.25 * 24 * 60 * 60
    if time_units == "years"
        times_sec = times .* yr2sec
    else
        times_sec = times
    end

    # Create figure
    fig = Figure(resolution=figsize)

    ax = Axis(fig[1, 1],
             xlabel="Accumulated Slip (m)",
             ylabel="Depth (km)",
             yreversed=true,
             title="Cumulative Slip Profiles at Specific Times")

    colors = [:blue, :red, :green, :purple, :orange, :cyan, :magenta, :brown]

    # Plot slip profile for each requested time
    for (i, t_requested) in enumerate(times_sec)
        # Find closest time index
        idx = argmin(abs.(times_global .- t_requested))
        actual_time = times_global[idx]

        # Read slip at all depths for this time step
        slip_profile = zeros(length(available_depths))

        for (j, depth) in enumerate(available_depths)
            time, slip, _, _, _ = read_depth_timeseries(h5_file, depth)
            slip_profile[j] = slip[idx]
        end

        # Filter by max depth
        if !isnothing(max_depth)
            mask = available_depths .<= max_depth
            slip_profile = slip_profile[mask]
            depths_plot = available_depths[mask]
        else
            depths_plot = available_depths
        end

        color = colors[mod1(i, length(colors))]
        label = if time_units == "years"
            @sprintf("%.1f yr", actual_time / yr2sec)
        else
            @sprintf("%.2e s", actual_time)
        end

        lines!(ax, slip_profile, depths_plot, color=color, linewidth=2, label=label)
    end

    # Add legend
    axislegend(ax, position=:rb)

    # Save figure
    save(output_file, fig)

    @info "Slip at times plot saved" file=output_file n_times=length(times)

    return fig
end


"""
    plot_slip_evolution(h5_file::String, depth::Float64;
                       output_file="slip_evolution.png",
                       time_units="years",
                       figsize=(800, 400))

Plot cumulative slip evolution at a specific depth.

# Arguments
- `h5_file::String`: Path to HDF5 output file
- `depth::Float64`: Depth in km
- `output_file::String`: Output filename (default: "slip_evolution.png")
- `time_units::String`: "years" or "seconds" (default: "years")
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (800, 400))

# Example
```julia
plot_slip_evolution("outputs/strike_slip_benchmark.h5", 10.0)
```
"""
function plot_slip_evolution(h5_file::String, depth::Float64;
                            output_file="slip_evolution.png",
                            time_units="years",
                            figsize=(800, 400))
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
             ylabel="Cumulative Slip (m)",
             title="Cumulative Slip at $(depth) km depth")

    lines!(ax, times, slip, color=:blue, linewidth=2)

    # Save figure
    save(output_file, fig)

    @info "Slip evolution plot saved" file=output_file depth=depth

    return fig
end


"""
    plot_slip_deficit(h5_file::String, plate_velocity::Float64;
                     output_file="slip_deficit.png",
                     time_units="years",
                     figsize=(800, 400),
                     depth=10.0)

Plot slip deficit (accumulated plate motion - cumulative slip) at a specific depth.

# Arguments
- `h5_file::String`: Path to HDF5 output file
- `plate_velocity::Float64`: Plate velocity in m/s
- `output_file::String`: Output filename (default: "slip_deficit.png")
- `time_units::String`: "years" or "seconds" (default: "years")
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (800, 400))
- `depth::Float64`: Depth in km (default: 10.0)

# Example
```julia
plot_slip_deficit("outputs/strike_slip_benchmark.h5", 1e-9, depth=10.0)
```
"""
function plot_slip_deficit(h5_file::String, plate_velocity::Float64;
                          output_file="slip_deficit.png",
                          time_units="years",
                          figsize=(800, 400),
                          depth=10.0)
    # Read data for specific depth
    times, slip, slip_rate, stress, state = read_depth_timeseries(h5_file, depth)

    # Calculate accumulated plate motion
    plate_motion = plate_velocity .* times

    # Calculate slip deficit
    slip_deficit = plate_motion .- slip

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
             ylabel="Slip Deficit (m)",
             title="Slip Deficit at $(depth) km depth")

    lines!(ax, times, slip_deficit, color=:red, linewidth=2)
    hlines!(ax, [0.0], color=:black, linestyle=:dash, linewidth=1)

    # Save figure
    save(output_file, fig)

    @info "Slip deficit plot saved" file=output_file depth=depth

    return fig
end
