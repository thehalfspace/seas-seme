# Scripts Directory

Utility scripts for SEAS-SEME simulation workflows.

## Available Scripts

### `generate_mesh.jl`

Generate unstructured quadrilateral meshes for earthquake cycle simulations.

**Usage with configuration file (recommended):**
```bash
julia --project=. scripts/generate_mesh.jl config/mesh_2d_homogeneous.toml simulation_name
```

**Usage with command line parameters:**
```bash
julia --project=. scripts/generate_mesh.jl simulation_name --Lx 60e3 --refine 0.1e3
```

**Examples:**

Standard mesh using default config:
```bash
julia --project=. scripts/generate_mesh.jl config/mesh_2d_homogeneous.toml strike_slip_benchmark
```

Quick test mesh with finer resolution:
```bash
julia --project=. scripts/generate_mesh.jl test_mesh --refine 0.1e3
```

Large domain:
```bash
julia --project=. scripts/generate_mesh.jl large_domain --Lx 80e3 --Ly 80e3
```

**Output:**
- Mesh files are saved to `data/mesh/{simulation_name}/`
- Generated files: `unstructured.mesh`, `unstructured.control`, `unstructured.tec`
- Visualize mesh in `notebooks/01_mesh_visualization.jl`
- Reference the mesh file in your simulation config.toml

**Mesh Configurations:**
- Default config: `config/mesh_2d_homogeneous.toml`
- Create custom configs for specific geometries (e.g., narrow fault zones)

**See also:**
- Full documentation: Run `julia --project=. scripts/generate_mesh.jl --help`
- Mesh module: `src/Mesh/Unstructured.jl`

---

### `run_simulation.jl`

Run earthquake cycle simulations from configuration files.

**Usage:**
```bash
julia --project=. scripts/run_simulation.jl config/strike_slip_2d.toml
```

**Prerequisites:**
- Generate a mesh first using `generate_mesh.jl`
- Ensure the mesh file path in your config matches the generated mesh

**Example workflow:**

1. Generate mesh:
```bash
julia --project=. scripts/generate_mesh.jl config/mesh_2d_homogeneous.toml strike_slip_benchmark
```

2. Update simulation config to reference the mesh:
```toml
[mesh]
file = "data/mesh/strike_slip_benchmark/unstructured.mesh"
polynomial_degree = 4
```

3. Run simulation:
```bash
julia --project=. scripts/run_simulation.jl config/strike_slip_2d.toml
```

**Output:**
- Simulation results saved to `data/{simulation_name}/`
- HDF5 output: `data/{simulation_name}/outputs/{simulation_name}.h5`
- Log file: `data/{simulation_name}/output.log`
- Parameters: `data/{simulation_name}/params/`

**See also:**
- Visualization: `notebooks/02_parameters_setup.jl`, `notebooks/03_results_output.jl`
- Configuration docs: `src/Config/Parameters.jl`
