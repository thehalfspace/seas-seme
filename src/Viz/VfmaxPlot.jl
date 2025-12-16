"""
    VfmaxPlot

Plotting functions for maximum slip rate (Vfmax) time series.
"""

using CairoMakie
using Printf


"""
    plot_vfmax(h5_file::String;
              output_file="vfmax.png",
              time_units="years",
              figsize=(800, 400),
              log_scale=true,
              time_range=nothing,
              vfmax_range=nothing)

Plot maximum fault slip rate (Vfmax) vs time.

# Arguments
- `h5_file::String`: Path to HDF5 output file
- `output_file::String`: Output filename (default: "vfmax.png")
- `time_units::String`: Time units: "years" or "seconds" (default: "years")
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (800, 400))
- `log_scale::Bool`: Use log scale for y-axis (default: true)
- `time_range::Union{Tuple{Float64,Float64},Nothing}`: Time range to plot (default: nothing = all)
- `vfmax_range::Union{Tuple{Float64,Float64},Nothing}`: Vfmax range for y-axis (default: nothing = auto)

# Example
```julia
plot_vfmax("outputs/strike_slip_benchmark.h5")
plot_vfmax("outputs/strike_slip_benchmark.h5", time_units="seconds", log_scale=false)
plot_vfmax("outputs/strike_slip_benchmark.h5", time_range=(0.0, 100.0))
```
"""
function plot_vfmax(h5_file::String;
                   output_file="vfmax.png",
                   time_units="years",
                   figsize=(800, 400),
                   log_scale=true,
                   time_range=nothing,
                   vfmax_range=nothing)
    # Read data
    times, Vfmax = read_timeseries_data(h5_file)

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

    # Filter by time range if specified
    if !isnothing(time_range)
        mask = (times .>= time_range[1]) .& (times .<= time_range[2])
        times = times[mask]
        Vfmax = Vfmax[mask]
    end

    # Create figure
    fig = Figure(resolution=figsize)

    ax = Axis(fig[1, 1],
             xlabel=time_label,
             ylabel="Max. Slip Rate (m/s)",
             title="Maximum Fault Slip Rate")

    # Set y-scale
    if log_scale
        ax.yscale = log10
    end

    # Plot
    lines!(ax, times, Vfmax, color=:blue, linewidth=2)

    # Set y-limits if specified
    if !isnothing(vfmax_range)
        ylims!(ax, vfmax_range...)
    end

    # Save figure
    save(output_file, fig)

    @info "Vfmax plot saved" file=output_file

    return fig
end


"""
    plot_vfmax_comparison(h5_files::Vector{String}, labels::Vector{String};
                         output_file="vfmax_comparison.png",
                         time_units="years",
                         figsize=(800, 400),
                         log_scale=true,
                         time_range=nothing)

Plot and compare maximum slip rate from multiple simulations.

# Arguments
- `h5_files::Vector{String}`: Paths to HDF5 output files
- `labels::Vector{String}`: Labels for each simulation
- `output_file::String`: Output filename (default: "vfmax_comparison.png")
- `time_units::String`: Time units: "years" or "seconds" (default: "years")
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (800, 400))
- `log_scale::Bool`: Use log scale for y-axis (default: true)
- `time_range::Union{Tuple{Float64,Float64},Nothing}`: Time range to plot (default: nothing = all)

# Example
```julia
files = ["sim1.h5", "sim2.h5", "sim3.h5"]
labels = ["Simulation 1", "Simulation 2", "Simulation 3"]
plot_vfmax_comparison(files, labels)
```
"""
function plot_vfmax_comparison(h5_files::Vector{String}, labels::Vector{String};
                              output_file="vfmax_comparison.png",
                              time_units="years",
                              figsize=(800, 400),
                              log_scale=true,
                              time_range=nothing)
    if length(h5_files) != length(labels)
        error("Number of files and labels must match")
    end

    # Create figure
    fig = Figure(resolution=figsize)

    # Convert time units
    yr2sec = 365.25 * 24 * 60 * 60
    if time_units == "years"
        time_label = "Time (years)"
        time_conversion = 1.0 / yr2sec
    elseif time_units == "seconds"
        time_label = "Time (s)"
        time_conversion = 1.0
    else
        error("time_units must be 'years' or 'seconds'")
    end

    ax = Axis(fig[1, 1],
             xlabel=time_label,
             ylabel="Max. Slip Rate (m/s)",
             title="Maximum Fault Slip Rate Comparison")

    # Set y-scale
    if log_scale
        ax.yscale = log10
    end

    # Plot each simulation
    colors = [:blue, :red, :green, :purple, :orange, :cyan, :magenta, :brown]
    for (i, (file, label)) in enumerate(zip(h5_files, labels))
        times, Vfmax = read_timeseries_data(file)
        times = times .* time_conversion

        # Filter by time range if specified
        if !isnothing(time_range)
            mask = (times .>= time_range[1]) .& (times .<= time_range[2])
            times = times[mask]
            Vfmax = Vfmax[mask]
        end

        color = colors[mod1(i, length(colors))]
        lines!(ax, times, Vfmax, label=label, color=color, linewidth=2)
    end

    # Add legend
    axislegend(ax, position=:rt)

    # Save figure
    save(output_file, fig)

    @info "Vfmax comparison plot saved" file=output_file

    return fig
end


"""
    plot_vfmax_with_events(h5_file::String;
                          output_file="vfmax_events.png",
                          time_units="years",
                          figsize=(800, 400),
                          event_threshold=1e-3,
                          time_range=nothing)

Plot maximum slip rate with earthquake events highlighted.

Events are identified when Vfmax exceeds the event_threshold.

# Arguments
- `h5_file::String`: Path to HDF5 output file
- `output_file::String`: Output filename (default: "vfmax_events.png")
- `time_units::String`: Time units: "years" or "seconds" (default: "years")
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (800, 400))
- `event_threshold::Float64`: Threshold for identifying events in m/s (default: 1e-3)
- `time_range::Union{Tuple{Float64,Float64},Nothing}`: Time range to plot (default: nothing = all)

# Example
```julia
plot_vfmax_with_events("outputs/strike_slip_benchmark.h5", event_threshold=1e-2)
```
"""
function plot_vfmax_with_events(h5_file::String;
                               output_file="vfmax_events.png",
                               time_units="years",
                               figsize=(800, 400),
                               event_threshold=1e-3,
                               time_range=nothing)
    # Read data
    times, Vfmax = read_timeseries_data(h5_file)

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

    # Filter by time range if specified
    if !isnothing(time_range)
        mask = (times .>= time_range[1]) .& (times .<= time_range[2])
        times = times[mask]
        Vfmax = Vfmax[mask]
    end

    # Identify events
    event_mask = Vfmax .>= event_threshold
    event_times = times[event_mask]
    event_vfmax = Vfmax[event_mask]

    # Create figure
    fig = Figure(resolution=figsize)

    ax = Axis(fig[1, 1],
             xlabel=time_label,
             ylabel="Max. Slip Rate (m/s)",
             yscale=log10,
             title="Maximum Fault Slip Rate with Events")

    # Plot slip rate
    lines!(ax, times, Vfmax, color=:blue, linewidth=2, label="Vfmax")

    # Highlight events
    if length(event_times) > 0
        scatter!(ax, event_times, event_vfmax, color=:red, markersize=8, label="Events")
    end

    # Add threshold line
    hlines!(ax, [event_threshold], color=:red, linestyle=:dash, linewidth=1,
           label="Event Threshold")

    # Add legend
    axislegend(ax, position=:rt)

    # Save figure
    save(output_file, fig)

    @info "Vfmax with events plot saved" file=output_file n_events=length(event_times)

    return fig
end
