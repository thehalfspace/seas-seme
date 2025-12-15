"""
    Simulation

Main simulation orchestration - builds and configures all components.

Provides the `build_simulation` function that constructs a complete
simulation from a configuration file, ready to run.
"""

"""
    Simulation{T}

Container for all simulation components.

# Fields
- `config::SimulationConfig`: Configuration parameters
- `mesh::UnstructuredSEMesh{T}`: Spectral element mesh
- `physics::MaterialProperties{T}`: Material properties
- `ics`: Initial conditions
- `params`: Simulation parameters (Vpl, V₀, f₀, etc.)
- `state::SimulationState{T}`: Current simulation state
- `qs_solver::QuasistaticSolver`: Quasi-static solver with AMG
- `dyn_solver::DynamicSolver{T}`: Dynamic solver
- `timestepper::AdaptiveTimestepper{T}`: Timestep controller
- `io_manager`: I/O manager for output
- `M_global::Vector{T}`: Global mass matrix (diagonal)
- `K_el::Array{T,3}`: Elemental stiffness matrices

# Usage
```julia
config = load_config("config.toml")
sim = build_simulation(config)
run!(sim)
```
"""
struct Simulation{T<:AbstractFloat}
    config::SimulationConfig
    mesh::UnstructuredSEMesh{T}
    physics::MaterialProperties{T}
    ics
    params
    state::SimulationState{T}
    qs_solver::QuasistaticSolver
    dyn_solver::DynamicSolver{T}
    timestepper::AdaptiveTimestepper{T}
    io_manager
    M_global::Vector{T}
    K_el::Array{T,3}
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
3. Compute mass and stiffness matrices
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
    println("\n" * "="^80)
    println("Building Simulation")
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

    # 3. Build mass and stiffness matrices
    println("\n[3/8] Computing elemental matrices...")
    basis = LobattoLegendreBasis(mesh.polynomial_degree)

    M_el, M_global = build_mass_matrices(mesh, physics, basis)
    K_el = build_stiffness_matrices(mesh, physics, basis)
    @printf("  Mass and stiffness matrices computed\n")

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

    ics = build_initial_conditions(config.physics, mesh)
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
        K_el, mesh.dof_id, mesh, fltni,
        tolerance=T(config.solvers.quasistatic.tolerance),
        max_iterations=config.solvers.quasistatic.max_iterations,
        amg_max_levels=config.solvers.quasistatic.amg_max_levels,
        verbose=true
    )
    @printf("  Quasi-static solver built (AMG-CG)\n")

    # Dynamic solver
    dyn_solver = DynamicSolver(dt_min, verbose=false)
    @printf("  Dynamic solver built (leap-frog)\n")

    # 7. Build timestepper
    println("\n[7/8] Configuring timestepper...")
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

    # 8. Initialize simulation state
    println("\n[8/8] Initializing simulation state...")
    state = SimulationState(mesh, ics, params, v_init=T(5.0e-4))
    @printf("  Simulation state initialized\n")
    @printf("  Initial max slip rate: %.3e m/s\n", maximum_fault_slip_rate(state))

    # 9. Save initial parameters
    println("\n[9/9] Saving initial parameters...")
    save_initial_parameters(config, ics, mesh)
    @printf("  Parameters saved to: %s\n", joinpath(config.simulation.output_dir, "params"))

    # 10. Create I/O manager
    io_manager = create_io_manager(config, mesh)

    println("\n" * "="^80)
    println("Simulation built successfully!")
    println("="^80)
    println("\nSimulation directory: $(config.simulation.output_dir)")
    println("├── config.toml")
    println("├── params/")
    println("│   ├── friction_parameters.dat")
    println("│   ├── initial_stress.dat")
    println("│   └── fault_coordinates.dat")
    println("├── outputs/")
    println("│   └── $(config.simulation.name).h5")
    println("└── checkpoints/")
    println("="^80)

    return Simulation{T}(
        config, mesh, physics, ics, params,
        state, qs_solver, dyn_solver, timestepper,
        io_manager, M_global, K_el
    )
end
