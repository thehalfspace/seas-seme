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

# ╔═╡ 9edf9ace-0662-4552-ba6c-5ef95ed3a6f6
begin
	#dat = "../data/plane_strain_2d/outputs/plane_strain_2d.h5"
	#out = "../data/strike_slip_2d/slip2.png"
	dat = "../data/strike_slip_2d/outputs/strike_slip_2d.h5"
	out = "../data/strike_slip_2d/slip2.png"
	plot_slip_contours(dat,
                 output_file=out,
                 max_depth=20.0,  # Optional: limit depth
				 dynamic_step=10,
				 #yreversed=true
	)
end

# ╔═╡ cdc9b7d6-5229-44fc-8cec-da4fd287db76
#=begin
	using HDF5, DelimitedFiles
	# Read all snapshot data
	times, slip, slip_rate, stress, state, depths_km = read_snapshots("../data/strike_slip_2d/outputs/strike_slip_2d.h5")

	config = get_snapshot_config("../data/strike_slip_2d/outputs/strike_slip_2d.h5")
end
=#

# ╔═╡ 7c4cacb0-28c3-42e9-ad9f-c88bb67c3485
#=begin
	plot_cumulative_slip("../data/strike_slip_2d_highres/outputs/strike_slip_2d_highres.h5",
                 output_file="../data/strike_slip_2d/slip.png",
				 dynamic_interval=0.1,
				 quasistatic_interval=2.0,
				 figsize=(800,600),
	)
end
=#

# ╔═╡ 8a6de99c-c0f5-41b5-883b-cc1cc6c8096d
#slip

# ╔═╡ f0618e41-c884-416a-8710-3d380fa0e2a2
#= begin
	using CairoMakie
	fig = Figure()
	ax = Axis(fig[1,1])
	idx1 = 1
	idx2 = 125
	for i in idx1:idx2 
		lines!(slip[:,i], linewidth=0.5)
	end
	fig
end
=#

# ╔═╡ 758ec0dd-9b0c-4f10-ab2f-74f137c4d7ed


# ╔═╡ Cell order:
# ╠═77811d42-e426-11f0-0cff-ade0daa50c92
# ╠═9edf9ace-0662-4552-ba6c-5ef95ed3a6f6
# ╠═cdc9b7d6-5229-44fc-8cec-da4fd287db76
# ╠═7c4cacb0-28c3-42e9-ad9f-c88bb67c3485
# ╠═8a6de99c-c0f5-41b5-883b-cc1cc6c8096d
# ╠═f0618e41-c884-416a-8710-3d380fa0e2a2
# ╠═758ec0dd-9b0c-4f10-ab2f-74f137c4d7ed
