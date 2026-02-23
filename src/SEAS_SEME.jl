"""
SEAS_SEME: Spectral Element Method for earthquake cycle Simulations

A modular, HPC-ready Julia package for simulating earthquake cycles using
the Spectral Element Method (SEM) with rate-state friction.

# Main modules:
- `Config`: TOML-based configuration management
- `Mesh`: Mesh generation and geometry
- `Physics`: Physical models (rate-state friction, initial conditions)
- `Discretization`: SEM discretization (mass, stiffness, matrix-free operators)
- `Solvers`: Time integration solvers (quasistatic, dynamic, fault boundary)
- `Simulation`: Main simulation orchestration
- `IO`: HDF5 output, checkpointing, logging
- `Viz`: Visualization tools for simulation results
"""
module SEAS_SEME

# Standard library imports
using LinearAlgebra
using SparseArrays
using Printf
using TOML
using Dates

# External dependencies
using Trixi
using HOHQMesh
using AlgebraicMultigrid
using IterativeSolvers
using LinearSolve
using StaticArrays
using StatsBase
using DelimitedFiles
using HDF5
using JLD2

# GPU dependencies (loaded conditionally; CUDA.jl will be a no-op if no GPU present)
using CUDA
using CUDA.CUSPARSE
using KernelAbstractions
using AMGX

# Include submodules
include("NumericalConstants.jl")

include("Config/Parameters.jl")
include("Config/Config.jl")

include("Mesh/Geometry.jl")
include("Mesh/Faces.jl")
include("Mesh/GeometricChecks.jl")
include("Mesh/Connectivity.jl")
include("Mesh/Boundaries.jl")
include("Mesh/Unstructured.jl")
include("Mesh/Mesh.jl")

include("Physics/RateStateFriction.jl")
include("Physics/InitialConditions.jl")
include("Physics/Physics.jl")

include("Discretization/MassMatrix.jl")
include("Discretization/StiffnessMatrix.jl")
include("Discretization/StiffnessMatrixPlaneStrain.jl")
include("Discretization/MassMatrixPlaneStrain.jl")
include("Discretization/MetricWeights.jl")
include("GPU/GPUKernels.jl")
include("Discretization/MatrixFreeOperator.jl")
include("Discretization/Discretization.jl")

include("Simulation/SimulationState.jl")
include("Simulation/SimulationStatePlaneStrain.jl")

include("Solvers/FaultSolvers.jl")
include("Solvers/FaultSolversPlaneStrain.jl")
include("Solvers/QuasistaticSolver.jl")
include("Solvers/QuasistaticSolverPlaneStrain.jl")
include("Solvers/DynamicSolver.jl")
include("Solvers/DynamicSolverPlaneStrain.jl")
include("Solvers/Timestepping.jl")
include("Solvers/Solvers.jl")

include("Simulation/TimeLoop.jl")
include("Simulation/Simulation.jl")

include("IO/HDF5Writer.jl")
include("IO/Checkpointing.jl")
include("IO/Logging.jl")
include("IO/IO.jl")

include("Viz/Viz.jl")

# Public API exports
export SimulationConfig, MeshConfig, PhysicsConfig, SolverConfig, OutputConfig, CheckpointConfig
export load_config, validate_config
export build_simulation, run!
export load_checkpoint, save_checkpoint
export AbstractSimulationState, SimulationState, SimulationStatePlaneStrain
export QuasistaticSolverGPU


# ============================================================================
# AMGX lifecycle: initialize once at module load, finalize at exit
# ============================================================================
function __init__()
    if CUDA.functional()
        try
            AMGX.initialize()
            AMGX.initialize_plugins()
            atexit(() -> begin
                try
                    AMGX.finalize_plugins()
                    AMGX.finalize()
                catch
                end
            end)
        catch e
            @warn "AMGX initialization failed; GPU QS solver unavailable" exception=e
        end
    end
end

end # module SEAS_SEME
