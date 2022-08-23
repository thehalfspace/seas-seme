### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 8e673bdb-0187-4c67-af6c-315acfdec147
begin
    import Pkg
    # careful: this is _not_ a reproducible environment
    # activate the global environment
    Pkg.activate()

    using CairoMakie, LinearAlgebra, Trixi, DelimitedFiles
end

# ╔═╡ dc89d92d-0569-44c1-95b6-12cad091ee39
# Main workflow
include("proc/time_loop.jl")

# ╔═╡ 13a26ae6-fd44-11ec-3ff6-e75b0d9328f1
# Step 1: Meshing (pre/mesh.jl)
# Set the directories in mesh.jl and time_loop.jl

# ╔═╡ 7cc7a591-2485-4fa3-b9d9-2158b4d3ae3d


# ╔═╡ 1bc3c5f4-3821-4411-9202-9a217e03ff7d
begin
	mesh_file = joinpath(mesh_path, "unstructured.m")
	mesh = UnstructuredMesh2D(mesh_file)
end

# ╔═╡ 30ed9afc-e004-4a60-9672-97d7e57562af
begin
	coords = mesh.corners
	connectivity = mesh.element_node_ids
	color = "blue"
end

# ╔═╡ b90e7fda-6938-42bc-aa34-6667d6c513d9
begin
	scene = CairoMakie.mesh(coords, connectivity, shading = false)
	wireframe!(scene[end][1], color = (:black, 0.6), linewidth = 3)
end

# ╔═╡ b0be7c81-2a9e-482d-8cff-c51600c5ae51
CairoMakie.mesh

# ╔═╡ Cell order:
# ╠═8e673bdb-0187-4c67-af6c-315acfdec147
# ╠═13a26ae6-fd44-11ec-3ff6-e75b0d9328f1
# ╠═dc89d92d-0569-44c1-95b6-12cad091ee39
# ╠═7cc7a591-2485-4fa3-b9d9-2158b4d3ae3d
# ╠═1bc3c5f4-3821-4411-9202-9a217e03ff7d
# ╠═30ed9afc-e004-4a60-9672-97d7e57562af
# ╠═b90e7fda-6938-42bc-aa34-6667d6c513d9
# ╠═b0be7c81-2a9e-482d-8cff-c51600c5ae51
