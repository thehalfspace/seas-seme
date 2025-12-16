# SEAS-SEME Visualization Notebooks

Interactive Pluto.jl notebooks for visualizing SEAS-SEME simulation data.

## Notebooks

### 01. Mesh Visualization (`01_mesh_visualization.jl`)

**Status**: Placeholder (to be implemented)

Will visualize:
- Spectral element mesh structure
- Fault boundary geometry
- Element connectivity
- Node distribution

### 02. Parameters and Setup (`02_parameters_setup.jl`)

Visualizes simulation parameters and initial conditions:

**Plots included:**
1. **Initial Conditions** - Dual-axis plot showing stresses (σ_n, τ) and friction parameter (a-b)
2. **Friction Parameters** - Three-panel plot of a, b, and Lc
3. **(a-b) Difference** - Showing velocity-weakening (VW) and velocity-strengthening (VS) regions
4. **Stress Distribution** - Normal and shear stress vs depth

**Features:**
- Interactive plotting with CairoMakie
- Summary statistics of all parameters
- Color-coded visualization of VW/VS regions
- Automatic plot saving to `plots/` directory

**Required data:**
- `data/{simulation_name}/params/friction_parameters.dat`
- `data/{simulation_name}/params/initial_stress.dat`

### 03. Results Output (`03_results_output.jl`)

Visualizes simulation results from HDF5 output files:

**Plots included:**
1. **Maximum Slip Rate (Vfmax)** - Time series showing earthquake cycles
2. **Vfmax with Event Detection** - Highlights earthquakes above threshold
3. **Earthquake Cycle Heatmap** - Slip rate vs depth and time
4. **Combined Vfmax + Heatmap** - Two-panel aligned view
5. **Slip Rate at Specific Depth** - Time series at chosen depth
6. **Cumulative Slip Evolution** - Cumulative slip vs time

**Features:**
- Interactive HDF5 data reading
- Log-scale visualizations for slip rates
- Customizable depth ranges
- Event detection and highlighting
- Summary statistics

**Required data:**
- `data/{simulation_name}/outputs/{simulation_name}.h5`

## Usage

### Running a Notebook

1. **Install Pluto.jl** (if not already installed):
   ```julia
   using Pkg
   Pkg.add("Pluto")
   ```

2. **Start Pluto**:
   ```julia
   using Pluto
   Pluto.run()
   ```

3. **Open a notebook** from the Pluto interface, or directly:
   ```julia
   Pluto.run(notebook="notebooks/02_parameters_setup.jl")
   ```

### Customizing Paths

Each notebook has a section to configure paths:

```julia
# Set simulation directory
simulation_name = "strike_slip_benchmark"  # Change this
base_dir = joinpath(dirname(pwd()), "data")
sim_dir = joinpath(base_dir, simulation_name)
```

Simply change `simulation_name` to match your simulation directory.

## Dependencies

All notebooks use:
- **CairoMakie.jl** - Publication-quality plotting
- **HDF5.jl** - Reading simulation output (03 only)
- **DelimitedFiles.jl** - Reading parameter files (02 only)
- **Printf.jl** - Formatting

These will be automatically installed by Pluto when you first run the notebooks.

## Output

Plots are automatically saved to:
```
plots/{simulation_name}/
├── initial_conditions.png
├── friction_parameters.png
├── ab_difference.png
├── stress_distribution.png
├── vfmax.png
├── vfmax_events.png
├── eq_cycle_heatmap.png
├── vfmax_with_heatmap.png
├── sliprate_depth_10.0km.png
└── cumulative_slip_10.0km.png
```

## Tips

1. **Interactive editing**: Pluto notebooks are reactive - changing a variable automatically updates all dependent cells
2. **Cell order**: Cells can be arranged in any order; Pluto automatically determines dependencies
3. **Export**: Use the export menu to save notebooks as HTML or PDF
4. **Live reload**: Notebooks automatically reload when the source file changes

## Troubleshooting

**Issue**: Plots not showing
- **Solution**: Make sure CairoMakie is properly installed and the theme is set

**Issue**: File not found errors
- **Solution**: Check that the paths in the configuration section match your directory structure

**Issue**: HDF5 read errors
- **Solution**: Verify the HDF5 file exists and is not corrupted. Check that the dataset paths match the expected structure.

## Next Steps

For more advanced visualization or scripting, use the `Viz` module directly:

```julia
using SEAS_SEME.Viz

plot_initial_conditions("data/strike_slip_benchmark/params")
plot_vfmax("data/strike_slip_benchmark/outputs/strike_slip_benchmark.h5")
```

See `src/Viz/README.md` for full documentation of the Viz module.
