"""
Solvers module: Time integration for earthquake cycle simulations.

This module provides:
- **FaultSolvers**: Rate-state friction boundary conditions (Newton-Raphson)
- **QuasistaticSolver**: AMG-preconditioned conjugate gradient for equilibrium
- **DynamicSolver**: Explicit leap-frog integration for wave propagation
- **Timestepping**: Adaptive timestep computation with solver mode switching

# Solver Modes

## Quasi-static Mode
Used during interseismic period and slow nucleation phases.
- Solves equilibrium: K*u = -f
- Matrix-free CG with AMG preconditioner
- Adaptive timestep based on nucleation length
- Rate-state friction enforced via direct calculation

## Dynamic Mode
Used during rapid earthquake slip (seismic events).
- Explicit leap-frog time integration
- Fixed timestep from CFL condition
- Rate-state friction enforced via Newton-Raphson search
- Absorbing boundary conditions for artificial boundaries

# Automatic Mode Switching

Solver automatically switches based on maximum fault slip velocity:
- QS → Dynamic: when Vf_max ≥ 5 mm/s (configurable)
- Dynamic → QS: when Vf_max < 2 mm/s (configurable)

Hysteresis prevents rapid mode switching.

# References
- Kaneko, Y., Lapusta, N., & Ampuero, J. P. (2008). Spectral element modeling of
  spontaneous earthquake rupture on rate and state faults: Effect of velocity-
  strengthening friction at shallow depths. JGR, 113(B9).

- Kaneko, Y., Ampuero, J. P., & Lapusta, N. (2011). Spectral-element simulations of
  long-term fault slip: Effect of low-rigidity layers on earthquake-cycle dynamics.
  JGR, 116(B10).

- Lapusta, N., Rice, J. R., Ben-Zion, Y., & Zheng, G. (2000). Elastodynamic analysis
  for slow tectonic loading with spontaneous rupture episodes on faults with
  rate- and state-dependent friction. JGR, 105(B10), 23765-23789.
"""

# Module exports
export fault_slip_rate, state_time_evolution, nr_search
export QuasistaticSolver, build_quasistatic_solver, quasistatic_step!
export QuasistaticSolverGPU, build_quasistatic_solver_gpu, build_quasistatic_solver_gpu_plane_strain
export DynamicSolver, dynamic_step!
export AdaptiveTimestepper, compute_timestep, determine_solver_mode

# Plane-strain solver exports
export QuasistaticSolverPlaneStrain, build_quasistatic_solver_plane_strain
export DynamicSolverPlaneStrain
export fault_traction_from_force_plane_strain
export compute_stick_traction_plane_strain
export apply_fault_traction_plane_strain!
export compute_fault_impedance_plane_strain
export set_fault_velocity_plane_strain!
export get_fault_tangential_velocity


# ============================================================================
# GPU storage helper: convert a CPU array to the same storage type as a target
# ============================================================================

"""
    _as_storage(target::AbstractVector, src::Vector) -> storage-matched array

Return `src` converted to the same storage type as `target`.
- If target is a Vector, returns src as-is.
- If target is a CuVector, wraps src in CuArray.
Used for upload-on-demand after CPU NR loops.
"""
_as_storage(::Vector, src::Vector) = src
_as_storage(::CuVector, src::Vector) = CuArray(src)
