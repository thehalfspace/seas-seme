"""
    Simulation

Main simulation orchestration - builds and configures all components.

Provides the `build_simulation` function that constructs a complete
simulation from a configuration file, ready to run.
"""

"""
    Simulation{T, S, QS, DS}

Container for all simulation components. Parameterized on:
- `T`: Float type
- `S`: State type (SimulationState{T} or SimulationStatePlaneStrain{T})
- `QS`: Quasi-static solver type
- `DS`: Dynamic solver type

# Fields
- `config::SimulationConfig`: Configuration parameters
- `mesh::UnstructuredSEMesh{T}`: Spectral element mesh
- `physics::MaterialProperties{T}`: Material properties
- `ics`: Initial conditions
- `params`: Simulation parameters (Vpl, V₀, f₀, etc.)
- `state::S`: Current simulation state
- `qs_solver::QS`: Quasi-static solver
- `dyn_solver::DS`: Dynamic solver
- `timestepper::AdaptiveTimestepper{T}`: Timestep controller
- `io_manager`: I/O manager for output
- `log_io::IO`: Log file IO stream
- `M_global::Vector{T}`: Global mass matrix (diagonal)
- `weights`: Compact metric weight arrays (MetricWeightsAntiplane or MetricWeightsPlaneStrain)
- `H::Matrix{T}`: Derivative matrix (shared across elements)
- `Ht::Matrix{T}`: H' (pre-transposed)

# Usage
```julia
config = load_config("config.toml")
sim = build_simulation(config)
run!(sim)
```
"""
struct Simulation{T<:AbstractFloat, S, QS, DS, W}
    config::SimulationConfig
    mesh::UnstructuredSEMesh{T}
    physics::MaterialProperties{T}
    ics
    params
    state::S
    qs_solver::QS
    dyn_solver::DS
    timestepper::AdaptiveTimestepper{T}
    io_manager
    log_io::IO
    M_global::Vector{T}
    weights::W          # MetricWeightsAntiplane{T} or MetricWeightsPlaneStrain{T}
    H::Matrix{T}
    Ht::Matrix{T}
end


"""
    _finish_simulation(config, mesh, physics, ics, params, state, qs_solver, dyn_solver,
                       dt_min, M_global, weights, H, Ht, log_io, label)

Shared tail for both builder variants: configure timestepper, save parameters,
create IO manager, and print summary. Called after all formulation-specific steps.
"""
function _finish_simulation(config, mesh, physics, ics, params, state, qs_solver,
                             dyn_solver, dt_min, M_global, weights, H, Ht, log_io, label)
    T = eltype(M_global)

    # 7. Build timestepper
    println("\n[7] Configuring timestepper...")
    hcell = minimum(diff(mesh.boundaries.fault.coords[2,:]))
    timestepper = AdaptiveTimestepper(
        dt_min,
        T(config.solvers.dt_max),
        hcell,
        physics.μ,
        T(config.solvers.dynamic.velocity_threshold_qs_to_dyn),
        T(config.solvers.dynamic.velocity_threshold_dyn_to_qs)
    )
    @printf("  Adaptive timestepping configured\n")
    @printf("  dt_min: %.3e s, dt_max: %.3e s\n", dt_min, config.solvers.dt_max)

    # 8. Save initial parameters
    println("\n[8] Saving initial parameters...")
    save_initial_parameters(config, ics, mesh)
    @printf("  Parameters saved to: %s\n", joinpath(config.simulation.output_dir, "params"))

    # 9. Create I/O manager
    io_manager = create_io_manager(config, mesh)

    println("\n" * "="^80)
    println("Simulation built successfully! ($label)")
    println("="^80)
    println("\nSimulation directory: $(config.simulation.output_dir)")
    println("├── config.toml")
    println("├── output.log")
    println("├── params/")
    println("│   ├── friction_parameters.dat")
    println("│   ├── initial_stress.dat")
    println("│   └── fault_coordinates.dat")
    println("├── outputs/")
    println("│   └── $(config.simulation.name).h5")
    println("└── checkpoints/")
    println("\nMesh: $(config.mesh.file)")
    println("  $(mesh.n_elements) elements, polynomial degree $(mesh.polynomial_degree)")
    println("="^80)

    return Simulation(
        config, mesh, physics, ics, params,
        state, qs_solver, dyn_solver, timestepper,
        io_manager, log_io, M_global, weights, H, Ht
    )
end


"""
    build_simulation(config::SimulationConfig; T=Float64) -> Simulation{T}

Build complete simulation from configuration.

# Arguments
- `config::SimulationConfig`: Configuration loaded from TOML
- `T::Type`: Floating point type (default: Float64)

# Returns
- `Simulation{T}`: Complete simulation ready to run

# Process
1. Load mesh
2. Build mesh connectivity, geometry, boundaries
3. Compute mass matrices and metric weights
4. Generate initial conditions
5. Build solvers (QS with AMG, dynamic)
6. Initialize simulation state
7. Create I/O manager
8. Create timestepper

# Example
```julia
config = load_config("examples/config/strike_slip_2d.toml")
sim = build_simulation(config)
run!(sim)
```
"""
function build_simulation(config::SimulationConfig; T::Type=Float64)
    if config.physics.formulation == :plane_strain
        return build_simulation_plane_strain(config, T=T)
    else
        return build_simulation_antiplane(config, T=T)
    end
end


"""
    build_simulation_antiplane(config; T=Float64) -> Simulation

Build antiplane (SH) simulation using tensor-product metric weights.
"""
function build_simulation_antiplane(config::SimulationConfig; T::Type=Float64)
    # Setup logging to both console and file
    _, log_io = setup_logging(config)

    println("\n" * "="^80)
    println("Building Simulation (Antiplane / SH)")
    println("="^80)

    # 1. Build mesh (needs physics config for boundary impedance)
    println("\n[1/8] Building mesh...")
    mesh = build_mesh(config.mesh, config.physics)
    @printf("  Mesh loaded: %d elements, %d DOFs\n", mesh.n_elements, mesh.ndof)
    @printf("  Polynomial degree: p=%d\n", mesh.polynomial_degree)

    # 2. Material properties
    println("\n[2/8] Setting material properties...")
    T = Float64  # Use Float64 for all computations
    physics = MaterialProperties(
        T(config.physics.density),
        T(config.physics.shear_velocity)
    )
    @printf("  Density: %.1f kg/m³\n", physics.ρ)
    @printf("  Shear velocity: %.1f m/s\n", physics.vs)
    @printf("  Shear modulus: %.2e Pa\n", physics.μ)

    # 3. Build mass matrices and metric weights
    println("\n[3/8] Computing elemental matrices...")
    basis = LobattoLegendreBasis(mesh.polynomial_degree)

    M_el, M_global = build_mass_matrices(mesh, physics, basis)
    weights = build_metric_weights_antiplane(mesh, physics, basis)

    # Derivative matrix (shared across all elements)
    H  = Matrix{T}(basis.derivative_matrix')   # D^T
    Ht = Matrix{T}(basis.derivative_matrix)    # D (transpose of H)
    @printf("  Mass matrices and antiplane metric weights computed\n")

    # 4. Adjust mass matrix for absorbing boundaries
    println("\n[4/8] Applying boundary conditions...")
    dt_min = T(config.solvers.dt_min_factor) * minimum(diff(mesh.boundaries.fault.coords[2,:])) / physics.vs
    absorbing_id = mesh.boundaries.absorbing.node_ids
    absorb_matrix = mesh.boundaries.absorbing.matrix

    M_global[absorbing_id] .+= (T(0.5) * dt_min) .* absorb_matrix
    @printf("  Minimum timestep (CFL): %.3e s\n", dt_min)
    @printf("  Absorbing boundaries: %d nodes\n", length(absorbing_id))
    @printf("  Fault nodes: %d\n", length(mesh.boundaries.fault.node_ids))
    @printf("  Creep boundary nodes: %d\n", length(mesh.boundaries.creep.node_ids))

    # 5. Generate initial conditions
    println("\n[5/8] Generating initial conditions...")
    params = (
        Vpl = T(config.physics.plate_velocity),
        Vo = T(config.physics.reference_slip_rate),
        fo = T(config.physics.reference_friction),
        yr2sec = T(365.25 * 24 * 3600)
    )

    ics = build_initial_conditions(config.physics.initial_conditions.file, mesh,
                                   config.physics, config.physics.dip_angle)
    @printf("  Initial conditions generated\n")
    @printf("  Fault depth range: %.1f - %.1f km\n",
           minimum(mesh.boundaries.fault.coords[2,:]) / 1000,
           maximum(mesh.boundaries.fault.coords[2,:]) / 1000)

    # 6. Build solvers
    println("\n[6/8] Building solvers...")

    # Non-fault/creep DOF indices (free DOFs)
    interface_id = vcat(mesh.boundaries.creep.node_ids,
                       mesh.boundaries.fault.node_ids)
    fltni = collect(1:mesh.ndof)
    deleteat!(fltni, sort(unique(interface_id)))

    # Quasi-static solver with AMG preconditioner
    qs_solver = build_quasistatic_solver(
        weights, H, Ht, mesh.dof_id, mesh, fltni,
        tolerance=T(config.solvers.quasistatic.tolerance),
        max_iterations=config.solvers.quasistatic.max_iterations,
        amg_max_levels=config.solvers.quasistatic.amg_max_levels,
        verbose=true
    )
    @printf("  Quasi-static solver built (AMG-CG)\n")

    # Dynamic solver
    dyn_solver = DynamicSolver(dt_min, verbose=false)
    @printf("  Dynamic solver built (leap-frog)\n")

    # 6b. Initialize simulation state
    println("\n[6b] Initializing simulation state...")
    # v_init = Vpl/2 so that Vf = 2*v_init = Vpl (steady state, consistent with τ⁰)
    v_init_val = T(params.Vpl) / 2
    state = SimulationState(mesh, ics, params, v_init=v_init_val)
    @printf("  Simulation state initialized (v_init = %.3e m/s)\n", v_init_val)
    @printf("  Initial max slip rate: %.3e m/s\n", maximum_fault_slip_rate(state))

    return _finish_simulation(config, mesh, physics, ics, params, state,
                              qs_solver, dyn_solver, dt_min, M_global, weights, H, Ht,
                              log_io, "Antiplane")
end


"""
    build_simulation_plane_strain(config; T=Float64) -> Simulation

Build plane-strain (P-SV) simulation with 2-component displacement.
Uses tensor-product metric weights for efficient matvec.
"""
function build_simulation_plane_strain(config::SimulationConfig; T::Type=Float64)
    _, log_io = setup_logging(config)

    println("\n" * "="^80)
    println("Building Simulation (Plane-Strain / P-SV)")
    println("="^80)
    @printf("  Dip angle: %.1f°\n", config.physics.dip_angle)
    @printf("  Poisson's ratio: %.3f\n", config.physics.poisson_ratio)

    # 1. Build mesh
    println("\n[1/9] Building mesh...")
    mesh = build_mesh(config.mesh, config.physics)
    @printf("  Mesh loaded: %d elements, %d spatial DOFs, %d total DOFs\n",
           mesh.n_elements, mesh.ndof, 2 * mesh.ndof)
    @printf("  Polynomial degree: p=%d\n", mesh.polynomial_degree)

    # 2. Material properties (3-arg constructor: ρ, vs, ν)
    println("\n[2/9] Setting material properties...")
    T = Float64
    physics = MaterialProperties(
        T(config.physics.density),
        T(config.physics.shear_velocity),
        T(config.physics.poisson_ratio)
    )
    @printf("  Density: %.1f kg/m³\n", physics.ρ)
    @printf("  Shear velocity: %.1f m/s\n", physics.vs)
    @printf("  P-wave velocity: %.1f m/s\n", physics.vp)
    @printf("  Shear modulus: %.2e Pa\n", physics.μ)
    @printf("  Lamé λ: %.2e Pa\n", physics.λ)
    @printf("  Poisson's ratio: %.3f\n", physics.ν)

    # 3. Build mass matrices and metric weights
    println("\n[3/9] Computing elemental matrices (plane-strain)...")
    basis = LobattoLegendreBasis(mesh.polynomial_degree)

    M_el, M_global = build_mass_matrices_plane_strain(mesh, physics, basis)
    weights = build_metric_weights_plane_strain(mesh, physics, basis)

    # Derivative matrix (shared across all elements)
    H  = Matrix{T}(basis.derivative_matrix')   # D^T
    Ht = Matrix{T}(basis.derivative_matrix)    # D (transpose of H)
    @printf("  Plane-strain mass matrices and metric weights computed\n")

    ndof = mesh.ndof

    # 4. CFL timestep (based on P-wave velocity, which is faster)
    println("\n[4/9] Applying boundary conditions...")
    hcell = minimum(diff(mesh.boundaries.fault.coords[2,:]))
    dt_min = T(config.solvers.dt_min_factor) * hcell / physics.vp
    @printf("  Minimum timestep (CFL, Vp-based): %.3e s\n", dt_min)

    # Absorbing boundary mass correction
    absorbing_id = mesh.boundaries.absorbing.node_ids
    absorb_matrix = mesh.boundaries.absorbing.matrix

    # Apply to both x and y components
    for i in eachindex(absorbing_id)
        nid = absorbing_id[i]
        M_global[nid]        += T(0.5) * dt_min * absorb_matrix[i]
        M_global[ndof + nid] += T(0.5) * dt_min * absorb_matrix[i]
    end

    @printf("  Absorbing boundaries: %d nodes\n", length(absorbing_id))
    @printf("  Fault nodes: %d\n", length(mesh.boundaries.fault.node_ids))
    @printf("  Creep boundary nodes: %d\n", length(mesh.boundaries.creep.node_ids))

    # 5. Generate initial conditions
    println("\n[5/9] Generating initial conditions...")

    # Compute effective plate velocity based on loading direction
    Vpl_raw = config.physics.plate_velocity
    Vpl_eff = if config.physics.loading_direction == :dip_slip
        Vpl_raw * cosd(config.physics.dip_angle)
    else  # :strike_slip
        Vpl_raw
    end
    @printf("  Loading direction: %s\n", config.physics.loading_direction)
    @printf("  Plate velocity (raw): %.3e m/s\n", Vpl_raw)
    @printf("  Plate velocity (effective tangential): %.3e m/s\n", Vpl_eff)

    params = (
        Vpl = T(Vpl_eff),
        Vo = T(config.physics.reference_slip_rate),
        fo = T(config.physics.reference_friction),
        yr2sec = T(365.25 * 24 * 3600)
    )

    ics = build_initial_conditions(config.physics.initial_conditions.file, mesh,
                                   config.physics, config.physics.dip_angle)
    @printf("  Initial conditions generated\n")
    @printf("  Fault depth range: %.1f - %.1f km\n",
           minimum(mesh.boundaries.fault.coords[2,:]) / 1000,
           maximum(mesh.boundaries.fault.coords[2,:]) / 1000)

    # 6. Build solvers
    println("\n[6/9] Building solvers...")

    # Free DOF indices in 2*ndof space
    # Both x and y components of fault/creep nodes are prescribed
    interface_id_spatial = sort(unique(vcat(mesh.boundaries.creep.node_ids,
                                            mesh.boundaries.fault.node_ids)))

    # In 2*ndof space: fault/creep nodes for both components
    interface_id_2n = vcat(interface_id_spatial, ndof .+ interface_id_spatial)
    fltni = collect(1:2*ndof)
    deleteat!(fltni, sort(interface_id_2n))

    qs_solver = build_quasistatic_solver_plane_strain(
        weights, H, Ht, mesh.dof_id, mesh, fltni,
        tolerance=T(config.solvers.quasistatic.tolerance),
        max_iterations=config.solvers.quasistatic.max_iterations,
        amg_max_levels=config.solvers.quasistatic.amg_max_levels,
        verbose=true
    )
    @printf("  Plane-strain quasi-static solver built (AMG-CG)\n")

    dyn_solver = DynamicSolverPlaneStrain(dt_min, verbose=false)
    @printf("  Plane-strain dynamic solver built (leap-frog)\n")

    # 6b. Initialize simulation state
    println("\n[6b] Initializing simulation state...")
    # v_init = Vpl_eff/2 so that Vf = 2*v_init = Vpl_eff (steady state, consistent with τ⁰)
    v_init_val = T(Vpl_eff) / 2
    state = SimulationStatePlaneStrain(mesh, ics, params, v_init=v_init_val)
    @printf("  Plane-strain simulation state initialized (v_init = %.3e m/s)\n", v_init_val)
    @printf("  Initial max slip rate: %.3e m/s\n", maximum_fault_slip_rate(state))
    @printf("  Formulation: plane-strain, dip angle: %.1f°, loading: %s\n",
            config.physics.dip_angle, config.physics.loading_direction)

    return _finish_simulation(config, mesh, physics, ics, params, state,
                              qs_solver, dyn_solver, dt_min, M_global, weights, H, Ht,
                              log_io, "Plane-Strain")
end
