### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 77811d42-e426-11f0-0cff-ade0daa50c92
begin
	using Pkg
	Pkg.activate("/Users/pthakur8/workInbox/2026/projects/sem/seas-seme")
	using SEAS_SEME.Viz
end

# ╔═╡ cdc9b7d6-5229-44fc-8cec-da4fd287db76
begin
	# Read all snapshot data
	times, slip, slip_rate, stress, state, depths_km = read_snapshots("../data/strike_slip_2d/outputs/strike_slip_2d.h5")

	config = get_snapshot_config("../data/strike_slip_2d/outputs/strike_slip_2d.h5")
end

# ╔═╡ Cell order:
# ╠═77811d42-e426-11f0-0cff-ade0daa50c92
# ╠═cdc9b7d6-5229-44fc-8cec-da4fd287db76
