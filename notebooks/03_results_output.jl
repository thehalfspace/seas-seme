### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 87c0ec82-9bca-4bc5-b478-918c9e494776
begin
	using Pkg
	Pkg.activate("/Users/pthakur8/workInbox/2026/projects/sem/seas-seme")
end

# ╔═╡ b879d036-9c6a-4685-94f4-ccd290693c91
begin
	using CairoMakie
	using HDF5
	using Printf
	using SEAS_SEME
end

# ╔═╡ 35dd971d-0f49-4880-a623-ccc5242a960f
let
	using SEAS_SEME.Viz
	
	# Read all snapshot data
	times, slip, slip_rate, stress, state, depths_km = read_snapshots("../data/strike_slip_2d/outputs/strike_slip_2d.h5")

	config = get_snapshot_config("../data/plane_strain_2d/outputs/plane_strain_2d.h5")
end

# ╔═╡ 3193a07c-03b3-11ed-38c7-15818fc5ee26
md"""
# 03. Results Output Visualization

Visualization of simulation results from HDF5 output files.

This notebook visualizes:
- Maximum slip rate (Vfmax) time series
- Earthquake cycle heatmaps (slip rate vs depth and time)
- Cumulative slip profiles
- Slip rate at specific depths
"""

# ╔═╡ b6211b4e-e958-4ca1-add1-8da399a61973
md"""
## 📁 Set Paths

Specify the simulation directory and HDF5 output file location.
"""

# ╔═╡ aee6645e-f699-4a7c-8a05-46ee0b9fe545
begin
	# Set simulation directory
	simulation_name = "plane_strain_2d"
	base_dir = joinpath(dirname(pwd()), "data")
	sim_dir = joinpath(base_dir, simulation_name)
	h5_file = joinpath(sim_dir, "outputs", "$(simulation_name).h5")

	# Create output directory for plots
	fig_dir = joinpath(dirname(pwd()), "plots", simulation_name)
	mkpath(fig_dir)

	@info "Paths configured" h5_file fig_dir
end

# ╔═╡ 1c03f79e-5e6f-4066-b719-e1c4b6cea3d6
md"""
## 📊 Read Simulation Data

Read time series and geometry from HDF5 file.
"""

# ╔═╡ 2c03f79e-5e6f-4066-b719-e1c4b6cea3d7
begin
	# Read global time series
	times, Vfmax = h5open(h5_file, "r") do file
		t = read(file, "timeseries/time")
		v = read(file, "timeseries/max_slip_rate")
		return t, v
	end

	# Read fault geometry
	fault_x, fault_y, depths_km = h5open(h5_file, "r") do file
		x = read(file, "mesh/fault_x")
		y = read(file, "mesh/fault_y")
		d = read(file, "mesh/fault_depths_km")
		return x, y, d
	end

	# Convert time to years
	yr2sec = 365.25 * 24 * 60 * 60
	times_years = times ./ yr2sec

	@info "Data loaded" n_timesteps=length(times) n_fault_nodes=length(depths_km)
end

# ╔═╡ 3c03f79e-5e6f-4066-b719-e1c4b6cea3d8
md"""
## 🎨 Global Figure Settings

Configure plotting style and theme.
"""

# ╔═╡ 4c03f79e-5e6f-4066-b719-e1c4b6cea3d9
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
## 📈 Plot 1: Maximum Slip Rate (Vfmax)

Time series of maximum fault slip rate showing earthquake cycles.
"""

# ╔═╡ 95dd971d-0f49-4880-a623-ccc5242a959b
let
	figsize = (900, 400)
	fig = Figure(resolution=figsize)

	ax = Axis(fig[1, 1],
			 xlabel="Time (years)",
			 ylabel="Max. Slip Rate (m/s)",
			 yscale=log10,
			 title="Maximum Fault Slip Rate")

	lines!(ax, times_years, Vfmax, color=:blue, linewidth=2.5)

	save(joinpath(fig_dir, "vfmax.png"), fig)
	fig
end

# ╔═╡ 25dd971d-0f49-4880-a623-ccc5242a960e
md"""
## 📈 Plot 4: Cumulative Slip Evolution

Cumulative slip vs time at a specific depth.
"""

# ╔═╡ fa9df7ba-e4b6-4d96-9f57-205c3966e234
md"""
## 📊 Summary Statistics

Key statistics about the simulation results.
"""

# ╔═╡ 0a9df7ba-e4b6-4d96-9f57-205c3966e235
begin
	# Calculate statistics
	n_events = count(Vfmax .>= 1e-3)
	max_vfmax = maximum(Vfmax)
	min_vfmax = minimum(Vfmax)
	total_time_years = maximum(times_years)

	md"""
	### Simulation Summary

	**Time Range:**
	- Total simulation time: $(round(total_time_years, digits=2)) years

	**Maximum Slip Rate (Vfmax):**
	- Maximum: $(max_vfmax) m/s
	- Minimum: $(min_vfmax) m/s
	- Number of events (Vfmax ≥ 1e-3 m/s): $(n_events)

	**Fault Geometry:**
	- Depth range: $(minimum(depths_km)) to $(maximum(depths_km)) km
	- Number of fault nodes: $(length(depths_km))

	**Time Steps:**
	- Total timesteps: $(length(times))
	"""
end

# ╔═╡ Cell order:
# ╟─3193a07c-03b3-11ed-38c7-15818fc5ee26
# ╠═87c0ec82-9bca-4bc5-b478-918c9e494776
# ╠═b879d036-9c6a-4685-94f4-ccd290693c91
# ╟─b6211b4e-e958-4ca1-add1-8da399a61973
# ╠═aee6645e-f699-4a7c-8a05-46ee0b9fe545
# ╟─1c03f79e-5e6f-4066-b719-e1c4b6cea3d6
# ╠═2c03f79e-5e6f-4066-b719-e1c4b6cea3d7
# ╟─3c03f79e-5e6f-4066-b719-e1c4b6cea3d8
# ╠═4c03f79e-5e6f-4066-b719-e1c4b6cea3d9
# ╟─85dd971d-0f49-4880-a623-ccc5242a959a
# ╠═95dd971d-0f49-4880-a623-ccc5242a959b
# ╟─25dd971d-0f49-4880-a623-ccc5242a960e
# ╠═35dd971d-0f49-4880-a623-ccc5242a960f
# ╟─fa9df7ba-e4b6-4d96-9f57-205c3966e234
# ╟─0a9df7ba-e4b6-4d96-9f57-205c3966e235
