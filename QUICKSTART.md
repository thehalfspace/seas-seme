# SEAS-SEME Quickstart Guide

Get up and running with earthquake cycle simulations in 5 minutes.

## Prerequisites

- Julia 1.10 or later
- 4+ GB RAM
- Linux/macOS (Windows may work but untested)

## Installation

### 1. Clone and Setup

```bash
cd /path/to/seas-seme
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

This will install all dependencies from `Project.toml`.

### 3. Verify Installation

```julia
julia --project=.
julia> using SEAS_SEME
julia> # Should load without errors
```

## Quick Test Run

### 1. Check Example Configuration

The example config is at `examples/config/strike_slip_2d.toml`. Key parameters:

```toml
[simulation]
total_time = 6.31152e8  # 20 years in seconds

[mesh]
file = "./data/structured.mesh"  # ← Make sure this file exists!

[solvers]
dt_max = 8.64e6  # 100 days
```

### 2. Verify Mesh File

**Important:** Check that your mesh file exists:
```bash
ls -lh data/structured.mesh
```

If missing, you need to generate it first. The mesh should be in Trixi/HOHQMesh format.

### 3. Run a Short Test (Recommended First)

For initial testing, modify the config to run just 1 year:

```bash
# Create a test config
cp examples/config/strike_slip_2d.toml examples/config/test_short.toml
```

Edit `test_short.toml`:
```toml
[simulation]
total_time = 3.15576e7  # 1 year instead of 20
```

Then run:
```bash
julia --project=. examples/scripts/run_simulation.jl examples/config/test_short.toml
```

### 4. Expected Output

You should see:
```
================================================================================
SEAS-SEME: Spectral Element Method for Earthquake Cycle Simulations
================================================================================
Configuration: examples/config/test_short.toml
Julia threads: 1
================================================================================

Loading configuration...
✓ Configuration loaded

Building simulation...

[1/8] Building mesh...
  Mesh loaded: X elements, Y DOFs
  Polynomial degree: p=4

[2/8] Setting material properties...
  Density: 2670.0 kg/m³
  Shear velocity: 3464.0 m/s
  Shear modulus: 3.20e+10 Pa

[3/8] Computing elemental matrices...
  Mass and stiffness matrices computed

[4/8] Applying boundary conditions...
  Minimum timestep (CFL): X.XXe-XX s
  Absorbing boundaries: XX nodes
  Fault nodes: XX

[5/8] Generating initial conditions...
  Initial conditions generated

[6/8] Building solvers...
Building AMG preconditioner...
  AMG levels: X
  Coarsest level size: XXX
  Quasi-static solver built (AMG-CG)
  Dynamic solver built (leap-frog)

[7/8] Configuring timestepper...
  Adaptive timestepping configured

[8/8] Initializing simulation state...
  Simulation state initialized

✓ Simulation built successfully!

Starting time integration...

[Xs] iter=1, t=X.XXe+XX s (0.0%), dt=X.XXe+XX s, Vf_max=X.XXe-XX m/s, mode=quasistatic
...
```

## Output Files

After running, check:

### 1. HDF5 Output

```bash
ls -lh data/*.h5
```

Should see: `data/strike_slip_benchmark.h5` (or your simulation name)

### 2. Checkpoints (if enabled)

```bash
ls -lh checkpoints/
```

Should see: `checkpoint_iter_XXXXX.jld2` files

## Examining Results

### Quick Python Script

```python
import h5py
import matplotlib.pyplot as plt

# Open HDF5 file
with h5py.File('data/strike_slip_benchmark.h5', 'r') as f:
    # List structure
    print("File structure:")
    def print_structure(name):
        print(name)
    f.visit(print_structure)

    # Read time series
    time = f['timeseries/time'][:]
    max_slip_rate = f['timeseries/max_slip_rate'][:]

    # Read specific depth (e.g., 10 km)
    slip_rate_10km = f['timeseries/depth_10.0km/slip_rate'][:]

    # Plot
    plt.figure(figsize=(12, 4))

    plt.subplot(1, 2, 1)
    plt.semilogy(time / (365.25*24*3600), max_slip_rate)
    plt.xlabel('Time (years)')
    plt.ylabel('Max Slip Rate (m/s)')
    plt.title('Maximum Fault Slip Rate')

    plt.subplot(1, 2, 2)
    plt.semilogy(time / (365.25*24*3600), slip_rate_10km)
    plt.xlabel('Time (years)')
    plt.ylabel('Slip Rate at 10 km (m/s)')
    plt.title('Slip Rate at 10 km Depth')

    plt.tight_layout()
    plt.savefig('results.png', dpi=150)
    print("Plot saved to results.png")
```

### Using Julia

```julia
using HDF5
using Plots

h5open("data/strike_slip_benchmark.h5", "r") do f
    time = read(f["timeseries/time"]) / (365.25*24*3600)  # Convert to years
    Vf_max = read(f["timeseries/max_slip_rate"])

    plot(time, Vf_max,
         yscale=:log10,
         xlabel="Time (years)",
         ylabel="Max Slip Rate (m/s)",
         title="Earthquake Cycle",
         legend=false)

    savefig("slip_rate.png")
end
```

## Restarting from Checkpoint

If your simulation was interrupted:

```bash
# Restart from latest checkpoint
julia --project=. examples/scripts/restart.jl --latest checkpoints/

# Or from specific checkpoint
julia --project=. examples/scripts/restart.jl checkpoints/checkpoint_iter_50000.jld2
```

## Common Issues

### 1. Mesh File Not Found

**Error:**
```
ERROR: Mesh file not found: ./data/structured.mesh
```

**Solution:**
- Check the mesh file path in your config
- Generate mesh using HOHQMesh or Trixi utilities
- Or update config to point to correct mesh location

### 2. Package Not Found

**Error:**
```
ERROR: ArgumentError: Package X not found in current path.
```

**Solution:**
```julia
using Pkg
Pkg.add("X")
```

### 3. Out of Memory

**Error:**
```
ERROR: OutOfMemoryError()
```

**Solution:**
- Reduce simulation time in config
- Use smaller mesh
- Reduce output frequency
- Run on machine with more RAM

### 4. CG Solver Not Converging

**Output:**
```
CG iterations: 100, converged: false
```

**Solutions:**
- Increase `max_iterations` in config
- Decrease `tolerance` (make it less strict)
- Check mesh quality
- Verify boundary conditions

### 5. NaN/Inf Values

**Error:**
```
ERROR: Non-finite values detected
```

**Solutions:**
- Reduce initial slip rate (`v_init` in SimulationState.jl)
- Reduce `dt_max` in config
- Increase `dt_min_factor` (more stable CFL)
- Check initial stress conditions

### 6. Very Slow Performance

**Symptoms:**
- Many iterations per second initially, then slows down
- Spending lots of time in quasi-static solver

**Solutions:**
- This is expected during slow interseismic periods
- Use larger `dt_max` for quasi-static phases
- Check AMG preconditioner quality (should converge in <20 iterations)
- Profile code to find bottlenecks

## Debugging Tips

### 1. Enable Verbose Output

Modify solver construction in `build_simulation()` to enable verbose mode:
```julia
qs_solver = build_quasistatic_solver(..., verbose=true)
dyn_solver = DynamicSolver(dt_min, verbose=true)
```

### 2. Check Intermediate State

Add debug output in `TimeLoop.jl`:
```julia
if state.iteration == 100
    @show state.u[1:10]
    @show state.v[1:10]
    @show state.Vf[1:10]
end
```

### 3. Save More Frequent Checkpoints

Edit config:
```toml
[checkpointing]
enabled = true
interval = 3.1557e6  # 0.1 years instead of 1 year
keep_last = 10
```

### 4. Reduce Problem Size

For faster iteration during debugging:
- Use smaller mesh (fewer elements)
- Shorter simulation time (1 year)
- Larger timesteps (within stability limits)

## Performance Tips

### 1. Use Multiple Threads

```bash
export JULIA_NUM_THREADS=16
julia --project=. examples/scripts/run_simulation.jl config.toml
```

**Note:** Current implementation may not be fully threaded yet. Check `@threads` usage.

### 2. Optimize Julia

```bash
# Precompile everything first
julia --project=. -e 'using SEAS_SEME'

# Run with optimizations
julia --project=. -O3 examples/scripts/run_simulation.jl config.toml
```

### 3. Monitor Resource Usage

```bash
# In another terminal
watch -n 1 'ps aux | grep julia'

# Or use htop
htop
```

## Next Steps

1. **Validate Results:** Compare with original code output
2. **Benchmark Performance:** Profile hot loops
3. **Test Edge Cases:** Different mesh sizes, parameters
4. **Add Tests:** Unit tests for key functions

## Getting Help

If you encounter issues:

1. Check error messages carefully
2. Verify mesh and config files
3. Try shorter simulation time first
4. Check available memory
5. Look for emergency checkpoints in `checkpoints/`

## Expected First Run

For a 1-year test simulation on a typical mesh:
- **Setup time:** 5-30 seconds (AMG preconditioner)
- **Run time:** Varies widely depending on mesh size and dynamics
  - Interseismic (QS): ~0.1-1 iterations/second
  - Coseismic (dynamic): ~100-1000 iterations/second
- **Memory usage:** 1-4 GB
- **Output file:** 1-100 MB (depends on output frequency)

For the full 20-year benchmark:
- **Expected time:** Hours to days depending on hardware
- **Earthquakes:** Should see periodic events (check slip rate spikes)
- **Output file:** Hundreds of MB to GB

Good luck! 🚀
