"""
    ParametersPlot

Plotting functions for simulation parameters (friction, stress, initial conditions).
"""

using CairoMakie
using Printf


"""
    plot_initial_conditions(params_dir::String;
                           output_file="initial_conditions.png",
                           figsize=(800, 600),
                           max_depth=nothing)

Plot friction parameters and initial stresses vs depth.

Creates a dual-axis plot:
- Left axis: Normal stress σ_n and shear stress τ (MPa) vs depth
- Right axis: (a-b) friction parameter vs depth

# Arguments
- `params_dir::String`: Path to params directory
- `output_file::String`: Output filename (default: "initial_conditions.png")
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (800, 600))
- `max_depth::Union{Float64,Nothing}`: Maximum depth to plot in km (default: nothing = all)

# Example
```julia
plot_initial_conditions("data/strike_slip_benchmark/params")
plot_initial_conditions("data/strike_slip_benchmark/params", max_depth=20.0)
```
"""
function plot_initial_conditions(params_dir::String;
                                output_file="initial_conditions.png",
                                figsize=(800, 600),
                                max_depth=nothing)
    # Read data
    depths, a, b, Lc, σ_n, τ = read_initial_conditions(params_dir)

    # Filter by max depth if specified
    if !isnothing(max_depth)
        mask = depths .<= max_depth
        depths = depths[mask]
        a = a[mask]
        b = b[mask]
        σ_n = σ_n[mask]
        τ = τ[mask]
    end

    # Create figure with dual axes
    fig = Figure(resolution=figsize)

    # Left axis: Stresses
    ax1 = Axis(fig[1, 1],
              xlabel="Stresses (MPa)",
              ylabel="Depth (km)",
              yreversed=true)

    lines!(ax1, σ_n, depths, label="Normal Stress σₙ", color=:black, linewidth=2)
    lines!(ax1, τ, depths, label="Shear Stress τ", color=:black, linewidth=2, linestyle=:dash)

    axislegend(ax1, position=:rb)

    # Right axis: (a-b)
    ax2 = Axis(fig[1, 1],
              xlabel="(a-b)",
              xlabelcolor=:blue,
              xticklabelcolor=:blue,
              yaxisposition=:right,
              yreversed=true)

    a_minus_b = a .- b
    lines!(ax2, a_minus_b, depths, color=:blue, linewidth=2)

    # Hide left y-axis for ax2
    hideydecorations!(ax2)

    # Save figure
    save(output_file, fig)

    @info "Initial conditions plot saved" file=output_file

    return fig
end


"""
    plot_friction_parameters(params_dir::String;
                            output_file="friction_parameters.png",
                            figsize=(1000, 400),
                            max_depth=nothing)

Plot friction parameters (a, b, Lc) vs depth in separate panels.

# Arguments
- `params_dir::String`: Path to params directory
- `output_file::String`: Output filename (default: "friction_parameters.png")
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (1000, 400))
- `max_depth::Union{Float64,Nothing}`: Maximum depth to plot in km (default: nothing = all)

# Example
```julia
plot_friction_parameters("data/strike_slip_benchmark/params")
```
"""
function plot_friction_parameters(params_dir::String;
                                 output_file="friction_parameters.png",
                                 figsize=(1000, 400),
                                 max_depth=nothing)
    # Read data
    depths, a, b, Lc, σ_n, τ = read_initial_conditions(params_dir)

    # Filter by max depth if specified
    if !isnothing(max_depth)
        mask = depths .<= max_depth
        depths = depths[mask]
        a = a[mask]
        b = b[mask]
        Lc = Lc[mask]
    end

    # Create figure with 3 panels
    fig = Figure(resolution=figsize)

    # Panel 1: a parameter
    ax1 = Axis(fig[1, 1],
              xlabel="a",
              ylabel="Depth (km)",
              yreversed=true,
              title="Parameter a")
    lines!(ax1, a, depths, color=:blue, linewidth=2)

    # Panel 2: b parameter
    ax2 = Axis(fig[1, 2],
              xlabel="b",
              yreversed=true,
              title="Parameter b")
    lines!(ax2, b, depths, color=:red, linewidth=2)
    hideydecorations!(ax2)

    # Panel 3: Lc parameter
    ax3 = Axis(fig[1, 3],
              xlabel="Lc (m)",
              yreversed=true,
              title="Characteristic Length")
    lines!(ax3, Lc, depths, color=:green, linewidth=2)
    hideydecorations!(ax3)

    # Link y-axes
    linkyaxes!(ax1, ax2, ax3)

    # Save figure
    save(output_file, fig)

    @info "Friction parameters plot saved" file=output_file

    return fig
end


"""
    plot_stress_distribution(params_dir::String;
                            output_file="stress_distribution.png",
                            figsize=(800, 600),
                            max_depth=nothing)

Plot normal and shear stress distribution vs depth.

# Arguments
- `params_dir::String`: Path to params directory
- `output_file::String`: Output filename (default: "stress_distribution.png")
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (800, 600))
- `max_depth::Union{Float64,Nothing}`: Maximum depth to plot in km (default: nothing = all)

# Example
```julia
plot_stress_distribution("data/strike_slip_benchmark/params")
```
"""
function plot_stress_distribution(params_dir::String;
                                 output_file="stress_distribution.png",
                                 figsize=(800, 600),
                                 max_depth=nothing)
    # Read data
    depths, a, b, Lc, σ_n, τ = read_initial_conditions(params_dir)

    # Filter by max depth if specified
    if !isnothing(max_depth)
        mask = depths .<= max_depth
        depths = depths[mask]
        σ_n = σ_n[mask]
        τ = τ[mask]
    end

    # Create figure
    fig = Figure(resolution=figsize)

    ax = Axis(fig[1, 1],
             xlabel="Stress (MPa)",
             ylabel="Depth (km)",
             yreversed=true,
             title="Stress Distribution")

    lines!(ax, σ_n, depths, label="Normal Stress σₙ", color=:blue, linewidth=2)
    lines!(ax, τ, depths, label="Shear Stress τ", color=:red, linewidth=2)

    axislegend(ax, position=:rb)

    # Save figure
    save(output_file, fig)

    @info "Stress distribution plot saved" file=output_file

    return fig
end


"""
    plot_ab_difference(params_dir::String;
                      output_file="ab_difference.png",
                      figsize=(800, 600),
                      max_depth=nothing)

Plot (a-b) difference to show velocity-weakening/strengthening regions.

Negative (a-b) indicates velocity-weakening (seismogenic).
Positive (a-b) indicates velocity-strengthening (stable).

# Arguments
- `params_dir::String`: Path to params directory
- `output_file::String`: Output filename (default: "ab_difference.png")
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (800, 600))
- `max_depth::Union{Float64,Nothing}`: Maximum depth to plot in km (default: nothing = all)

# Example
```julia
plot_ab_difference("data/strike_slip_benchmark/params")
```
"""
function plot_ab_difference(params_dir::String;
                           output_file="ab_difference.png",
                           figsize=(800, 600),
                           max_depth=nothing)
    # Read data
    depths, a, b, Lc, σ_n, τ = read_initial_conditions(params_dir)

    # Filter by max depth if specified
    if !isnothing(max_depth)
        mask = depths .<= max_depth
        depths = depths[mask]
        a = a[mask]
        b = b[mask]
    end

    a_minus_b = a .- b

    # Create figure
    fig = Figure(resolution=figsize)

    ax = Axis(fig[1, 1],
             xlabel="(a-b)",
             ylabel="Depth (km)",
             yreversed=true,
             title="Friction Parameter (a-b)")

    # Color by velocity-weakening (VW) vs velocity-strengthening (VS)
    vw_mask = a_minus_b .< 0
    vs_mask = a_minus_b .>= 0

    if any(vw_mask)
        lines!(ax, a_minus_b[vw_mask], depths[vw_mask],
              label="Velocity-Weakening (VW)", color=:red, linewidth=2)
    end

    if any(vs_mask)
        lines!(ax, a_minus_b[vs_mask], depths[vs_mask],
              label="Velocity-Strengthening (VS)", color=:blue, linewidth=2)
    end

    # Add vertical line at (a-b) = 0
    vlines!(ax, [0.0], color=:black, linestyle=:dash, linewidth=1)

    axislegend(ax, position=:rb)

    # Save figure
    save(output_file, fig)

    @info "(a-b) difference plot saved" file=output_file

    return fig
end
