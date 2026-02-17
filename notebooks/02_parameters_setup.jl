### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ b879d036-9c6a-4685-94f4-ccd290693c91
begin
	using Pkg
	Pkg.activate("/Users/pthakur8/workInbox/2026/projects/sem/seas-seme")
	using SEAS_SEME.Viz

	using CairoMakie
	using DelimitedFiles
	using Printf
end

# ╔═╡ 3193a07c-03b3-11ed-38c7-15818fc5ee26
md"""
# 02. Parameters and Setup Visualization

Visualization of simulation parameters and initial conditions.

This notebook visualizes:
- Friction parameters (a, b, Lc)
- Initial stress distribution (σ_n, τ)
- (a-b) to identify velocity-weakening/strengthening regions
"""

# ╔═╡ d35d3229-33a6-4f04-b270-ad1783e4fe4b
"a" * "b"

# ╔═╡ b6211b4e-e958-4ca1-add1-8da399a61973
md"""
## 📁 Set Paths

Specify the simulation directory and parameter files location.
"""

# ╔═╡ aee6645e-f699-4a7c-8a05-46ee0b9fe545
begin
	# Set simulation directory
	simulation_name = "dip_slip_2d"
	base_dir = joinpath(dirname(pwd()), "data")
	sim_dir = joinpath(base_dir, simulation_name)
	params_dir = joinpath(sim_dir, "params")

	data_file = joinpath(sim_dir, "outputs", simulation_name * ".h5")

	# Create output directory for plots
	fig_dir = joinpath(dirname(pwd()), "plots", simulation_name)
	mkpath(fig_dir)

	@info "Paths configured" params_dir fig_dir
end

# ╔═╡ 05bb5cbd-e6db-45be-95c8-9ec5d4da6dbb
md"""
## 📊 Read Parameter Data

Read friction parameters and initial stresses from the params directory.
"""

# ╔═╡ 1c03f79e-5e6f-4066-b719-e1c4b6cea3d6
begin
	# Read friction parameters
	friction_file = joinpath(params_dir, "friction_parameters.dat")
	friction_data = readdlm(friction_file, skipstart=2)

	depths_friction = friction_data[:, 1]
	a = friction_data[:, 2]
	b = friction_data[:, 3]
	Lc = friction_data[:, 4]

	# Read initial stress
	stress_file = joinpath(params_dir, "initial_stress.dat")
	stress_data = readdlm(stress_file, skipstart=2)

	depths_stress = stress_data[:, 1]
	σ_n = stress_data[:, 2]
	τ = stress_data[:, 3]

	@info "Data loaded" n_points=length(depths_friction)
end

# ╔═╡ 0c03f79e-5e6f-4066-b719-e1c4b6cea3d7
md"""
## 🎨 Global Figure Settings

Configure plotting style and theme.
"""

# ╔═╡ 2c03f79e-5e6f-4066-b719-e1c4b6cea3d8
begin
	# Global figure properties
	dpi = 150
	my_theme = Theme(
		fontsize=16,
		linewidth=2.5,
		Axis = (
			xlabelsize = 14,
			ylabelsize = 14,
			titlesize = 16
		)
	)
	set_theme!(my_theme)
end

# ╔═╡ 85dd971d-0f49-4880-a623-ccc5242a959a
md"""
## 📈 Plot 1: Initial Conditions (Stresses + Friction)

Dual-axis plot showing stresses on the left and (a-b) friction parameter on the right.
"""

# ╔═╡ 95dd971d-0f49-4880-a623-ccc5242a959b
let
	figsize = (800, 600)
	fig = Figure(resolution=figsize)

	# Left axis: Stresses
	ax1 = Axis(fig[1, 1],
			  xlabel="Stresses (MPa)",
			  ylabel="Depth (km)",
			  yreversed=true,
			  title="Initial Conditions")

	lines!(ax1, σ_n, depths_stress, label="Normal Stress σₙ", color=:black, linewidth=2.5)
	lines!(ax1, τ, depths_stress, label="Shear Stress τ", color=:black, linewidth=2.5, linestyle=:dash)

	axislegend(ax1, position=:rb)

	# Right axis: (a-b)
	ax2 = Axis(fig[1, 1],
			  xlabel="(a-b)",
			  xlabelcolor=:blue,
			  xticklabelcolor=:blue,
			  yaxisposition=:right,
			  yreversed=true)

	a_minus_b = a .- b
	lines!(ax2, a_minus_b, depths_friction, color=:blue, linewidth=2.5)

	hideydecorations!(ax2)

	save(joinpath(fig_dir, "initial_conditions.png"), fig)
	fig
end

# ╔═╡ 0135772b-4ba6-435f-b2fb-abfcb573fe06


# ╔═╡ 300c13cb-4f04-4a59-8d10-753971620ade
begin
	plot_slip_contours(data_file,
                 output_file=joinpath(fig_dir, "cumulative_slip.png"),
				 dynamic_step=10,
				 #yreversed=true
	)
end

# ╔═╡ dabb56ff-65dc-48a9-b061-39615f6ba0b1
begin
	# Read all snapshot data
	times, slip, slip_rate, stress, state, depths_km = read_snapshots("../data/dip_slip_2d/outputs/dip_slip_2d.h5")

	#config = get_snapshot_config("../data/strike_slip_2d/outputs/strike_slip_2d.h5")
end


# ╔═╡ d51cc2b3-8457-4ecd-aa2a-27f977285326
depths_km

# ╔═╡ b697b655-8aa5-42e9-aeb7-f81e20970676
begin
	figsize = (800, 600)
	fig = Figure(resolution=figsize)
	x1 = Axis(fig[1, 1],
			  xlabel="stress",
			  ylabel="Depth (km)",
			  title="slip")
	plot!(stress')
	fig
end

# ╔═╡ fa9df7ba-e4b6-4d96-9f57-205c3966e234
md"""
## 📊 Summary Statistics

Key statistics about the simulation parameters.
"""

# ╔═╡ Cell order:
# ╟─3193a07c-03b3-11ed-38c7-15818fc5ee26
# ╠═b879d036-9c6a-4685-94f4-ccd290693c91
# ╠═d35d3229-33a6-4f04-b270-ad1783e4fe4b
# ╟─b6211b4e-e958-4ca1-add1-8da399a61973
# ╠═aee6645e-f699-4a7c-8a05-46ee0b9fe545
# ╟─05bb5cbd-e6db-45be-95c8-9ec5d4da6dbb
# ╠═1c03f79e-5e6f-4066-b719-e1c4b6cea3d6
# ╟─0c03f79e-5e6f-4066-b719-e1c4b6cea3d7
# ╟─2c03f79e-5e6f-4066-b719-e1c4b6cea3d8
# ╟─85dd971d-0f49-4880-a623-ccc5242a959a
# ╠═95dd971d-0f49-4880-a623-ccc5242a959b
# ╠═0135772b-4ba6-435f-b2fb-abfcb573fe06
# ╠═300c13cb-4f04-4a59-8d10-753971620ade
# ╠═dabb56ff-65dc-48a9-b061-39615f6ba0b1
# ╠═d51cc2b3-8457-4ecd-aa2a-27f977285326
# ╠═b697b655-8aa5-42e9-aeb7-f81e20970676
# ╟─fa9df7ba-e4b6-4d96-9f57-205c3966e234
