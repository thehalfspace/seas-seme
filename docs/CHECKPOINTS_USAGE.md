# Checkpointing and Simulation Resume

This document explains how checkpoints work and how to resume a simulation from one.

## Overview

SEAS-SEME saves periodic checkpoints using JLD2. Each checkpoint contains the full
simulation state (displacements, velocities, fault slip rates, state variables, solver
objects, and initial conditions), allowing exact resumption from any saved point.

Checkpoints are useful for:
- Resuming long simulations after crashes or time limits
- Testing fixes to the solver by restarting from a known failing state
- Branching from a common interseismic state to test different parameters

## Configuration

```toml
[checkpointing]
enabled = true
interval_years = 100.0   # Save every N simulated years
keep_last = 5            # Number of normal checkpoints to retain
```

Checkpoints are written to `<output_dir>/<simulation_name>/checkpoints/`.

Normal checkpoints: `checkpoint_iter_<N>.jld2`
Emergency checkpoints (saved on crash): `checkpoint_emergency_iter_<N>.jld2`

## Resuming from a Checkpoint

Use `scripts/resume_checkpoint.jl`:

```bash
julia --project=. scripts/resume_checkpoint.jl <config_file> <checkpoint_file>
```

The script rebuilds the simulation from the config (needed to reconstruct the mesh,
which is not fully serializable), then loads all state from the checkpoint and resumes
the time loop from the saved iteration.

### Example: Testing a post-dynamic fix

The dip-slip 2D simulation crashed at iter=187198, the first quasistatic step after
a dynamic→quasistatic mode switch. An emergency checkpoint was saved automatically.
To test the CG warm-start fix without re-running 5 hours of interseismic loading:

```bash
julia --project=. scripts/resume_checkpoint.jl \
    config/dip_slip_2d.toml \
    data/dip_slip_2d_test/checkpoints/checkpoint_emergency_iter_187198.jld2
```

This takes ~1-2 minutes to rebuild the AMG preconditioner, then immediately resumes
at the failing iteration to verify the fix.

## Notes

- Checkpoints are backward-incompatible when struct fields are added or reordered.
  Emergency checkpoints from a run using the current code can always be resumed with
  the same code version.
- The mesh is always rebuilt from the config file (not stored in the checkpoint) because
  it contains references to the HOHQMesh file that are cheaper to reconstruct.
- Emergency checkpoints are never rotated (all are kept regardless of `keep_last`).
