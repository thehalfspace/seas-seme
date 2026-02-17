#!/usr/bin/env julia
"""
Diagnostic script: Inspect fault_mat (boundary mass) distribution
and check for shared nodes between fault and creep boundaries.

Usage:
    julia scripts/diagnose_fault_mat.jl config/dip_slip_2d.toml
"""

using Pkg
Pkg.activate(".")

using SEAS_SEME

config_file = length(ARGS) >= 1 ? ARGS[1] : "config/dip_slip_2d.toml"

println("="^60)
println("Fault Boundary Mass Diagnostic")
println("="^60)
println("Config: $config_file")

config = SEAS_SEME.load_config(config_file)
mesh = SEAS_SEME.build_mesh(config.mesh, config.physics)

fault = mesh.boundaries.fault
creep = mesh.boundaries.creep

println("\n--- Boundary sizes ---")
println("  Fault nodes: $(length(fault.node_ids))")
println("  Creep nodes: $(length(creep.node_ids))")

# Check for shared nodes between fault and creep
shared = intersect(Set(fault.node_ids), Set(creep.node_ids))
println("\n--- Shared fault/creep nodes ---")
println("  Number of shared nodes: $(length(shared))")
if !isempty(shared)
    for nid in shared
        fi = findfirst(==(nid), fault.node_ids)
        ci = findfirst(==(nid), creep.node_ids)
        println("    Node $nid: fault_mat=$(fault.matrix[fi]), creep_mat=$(creep.matrix[ci])")
        println("      fault coords: ($(fault.coords[1,fi]), $(fault.coords[2,fi]))")
    end
end

# Check for shared nodes between fault and absorbing
absorb = mesh.boundaries.absorbing
shared_fa = intersect(Set(fault.node_ids), Set(absorb.node_ids))
println("\n--- Shared fault/absorbing nodes ---")
println("  Number of shared nodes: $(length(shared_fa))")

# Fault_mat statistics
fm = fault.matrix
println("\n--- Fault matrix (boundary mass) statistics ---")
println("  min:    $(minimum(fm))")
println("  max:    $(maximum(fm))")
println("  mean:   $(sum(fm)/length(fm))")
println("  median: $(sort(fm)[div(length(fm),2)])")
println("  ratio max/min: $(maximum(fm)/minimum(fm))")

# Sort and show distribution
sorted_idx = sortperm(fm)
println("\n--- Smallest 10 fault_mat values ---")
for k in 1:min(10, length(fm))
    i = sorted_idx[k]
    x, y = fault.coords[1,i], fault.coords[2,i]
    println("  [$k] node_id=$(fault.node_ids[i]), fault_mat=$(fm[i]), coords=($x, $y)")
end

println("\n--- Largest 10 fault_mat values ---")
for k in 0:min(9, length(fm)-1)
    i = sorted_idx[end-k]
    x, y = fault.coords[1,i], fault.coords[2,i]
    println("  [$(length(fm)-k)] node_id=$(fault.node_ids[i]), fault_mat=$(fm[i]), coords=($x, $y)")
end

# Check tangent/normal vectors at endpoints
println("\n--- Fault endpoint geometry ---")
for label in ["First", "Last"]
    i = label == "First" ? 1 : length(fault.node_ids)
    x, y = fault.coords[1,i], fault.coords[2,i]
    tx, ty = fault.tangent[1,i], fault.tangent[2,i]
    nx, ny = fault.normal[1,i], fault.normal[2,i]
    println("  $label node (i=$i): coords=($x, $y)")
    println("    tangent=($tx, $ty), normal=($nx, $ny)")
    println("    fault_mat=$(fm[i])")
end

# GLL weights for reference
using SEAS_SEME: LobattoLegendreBasis
p = config.mesh.polynomial_degree
basis = LobattoLegendreBasis(p)
println("\n--- GLL weights (p=$p) ---")
println("  weights = $(basis.weights)")
println("  w_corner = $(basis.weights[1]), w_max = $(maximum(basis.weights))")
println("  ratio w_max/w_corner = $(maximum(basis.weights)/basis.weights[1])")

# Expected fault_mat range for uniform elements
# For fault edge length h, J_1D = h/2
# Single element contribution: w_k * h/2
# Corner node (shared by 2 elements): 2 * w_corner * h/2 = w_corner * h
# Interior node (single element):     w_interior * h/2
h_fault = config.mesh.polynomial_degree > 0 ? 1000.0 : 2000.0  # approximate from config
println("\n--- Expected fault_mat for h_fault ≈ $(h_fault) m ---")
println("  Corner shared (2 elem): $(2 * basis.weights[1] * h_fault/2)")
println("  Corner unshared (1 elem): $(basis.weights[1] * h_fault/2)")
println("  Interior (1 elem): $(maximum(basis.weights) * h_fault/2)")

# Active fault mask
mask = mesh.active_fault_mask
n_active = count(mask)
n_excluded = count(!, mask)
println("\n--- Active fault mask ---")
println("  Active: $n_active / $(length(mask))")
println("  Excluded: $n_excluded")
if n_excluded > 0
    println("  Excluded nodes:")
    for i in eachindex(fault.node_ids)
        if !mask[i]
            x, y = fault.coords[1,i], fault.coords[2,i]
            println("    [$i] node_id=$(fault.node_ids[i]), fault_mat=$(fm[i]), coords=($x, $y)")
        end
    end
end

println("\n" * "="^60)
println("Done.")
