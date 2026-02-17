"""
Main mesh module: defines UnstructuredSEMesh type and constructor.

This module ties together connectivity, geometry, and boundary extraction
to create a complete spectral element mesh representation.
"""

"""
    BoundaryInfo{T}

Container for all boundary data (fault, creep, absorbing, free surface).

# Fields
- `fault::BoundaryData{T}`: Fault boundary nodes and impedance
- `creep::BoundaryData{T}`: Creep boundary nodes and impedance
- `absorbing::BoundaryData{T}`: Absorbing boundary nodes and impedance
"""
struct BoundaryInfo{T<:AbstractFloat}
    fault::BoundaryData{T}
    creep::BoundaryData{T}
    absorbing::BoundaryData{T}
end

"""
    UnstructuredSEMesh{T}

Complete spectral element mesh with geometry, connectivity, and boundaries.

# Fields
- `trixi_mesh::UnstructuredMesh2D`: Underlying Trixi.jl mesh
- `dof_id::Array{Int,3}`: Connectivity matrix [nnodes, nnodes, nelements]
- `node_coords::Array{T,4}`: Node coordinates [2, nnodes, nnodes, nelements]
- `jac_matrix::Array{T,5}`: Jacobian matrices [2, 2, nnodes, nnodes, nelements]
- `jac_det::Array{T,3}`: Jacobian determinants [nnodes, nnodes, nelements]
- `normal_dirs::Array{T,4}`: Normal vectors [2, nnodes, 4, nelements]
- `boundaries::BoundaryInfo{T}`: Boundary node information
- `face_map::FaceMap`: Face abstraction (single source of truth for faces)
- `active_fault_mask::BitVector`: Which fault nodes participate in friction (true=active)
- `ndof::Int`: Total number of global DOFs
- `n_elements::Int`: Number of elements
- `polynomial_degree::Int`: Polynomial degree (p)
"""
struct UnstructuredSEMesh{T<:AbstractFloat}
    trixi_mesh::UnstructuredMesh2D
    dof_id::Array{Int,3}
    node_coords::Array{T,4}
    jac_matrix::Array{T,5}
    jac_det::Array{T,3}
    normal_dirs::Array{T,4}
    boundaries::BoundaryInfo{T}
    face_map::FaceMap
    active_fault_mask::BitVector
    ndof::Int
    n_elements::Int
    polynomial_degree::Int
end

"""
    build_mesh(config::MeshConfig, physics::PhysicsConfig) -> UnstructuredSEMesh

Construct a complete spectral element mesh from configuration.

# Arguments
- `config::MeshConfig`: Mesh configuration (file path, polynomial degree)
- `physics::PhysicsConfig`: Physics configuration (density, shear velocity for impedance)

# Returns
- `UnstructuredSEMesh{Float64}`: Complete mesh with geometry and boundaries

# Steps
1. Load HOHQMesh file using Trixi.jl
2. Build connectivity matrix (global DOF numbering)
3. Compute node coordinates in physical space
4. Compute Jacobian matrices and determinants
5. Compute normal vectors on boundaries
6. Extract boundary nodes (fault, creep, absorbing)

# Example
```julia
config = load_config("config.toml")
mesh = build_mesh(config.mesh, config.physics)
```
"""
function build_mesh(config::MeshConfig, physics::PhysicsConfig)::UnstructuredSEMesh{Float64}
    # Load mesh file
    @info "Loading mesh from: $(config.file)"
    if !isfile(config.file)
        error("Mesh file not found: $(config.file)")
    end

    trixi_mesh = UnstructuredMesh2D(config.file)

    # Check polynomial degree matches
    if trixi_mesh.polydeg != config.polynomial_degree
        @warn "Mesh polynomial degree ($(trixi_mesh.polydeg)) differs from config ($(config.polynomial_degree)). Using config value."
        # Note: May need to regenerate mesh with correct polynomial degree
    end

    polydeg = config.polynomial_degree
    nnodes = polydeg + 1
    n_elements = trixi_mesh.n_elements

    @info "Mesh: $(n_elements) elements, polynomial degree = $(polydeg)"

    # Get Gauss-Lobatto basis
    basis = LobattoLegendreBasis(polydeg)
    nodes = basis.nodes

    # 1. Build connectivity matrix
    @info "Building connectivity matrix..."
    dof_id = connectivity_matrix(trixi_mesh)
    ndof = maximum(dof_id)
    @info "Total DOFs: $(ndof)"

    # 2. Allocate geometry arrays
    node_coords = zeros(2, nnodes, nnodes, n_elements)
    jac_matrix = zeros(2, 2, nnodes, nnodes, n_elements)
    jac_det = zeros(nnodes, nnodes, n_elements)
    normal_dirs = zeros(2, nnodes, 4, n_elements)

    # 3. Compute geometry for each element
    @info "Computing geometry (coordinates, Jacobians, normals)..."
    for el in 1:n_elements
        # Get element corner points
        corners = trixi_mesh.corners[:, trixi_mesh.element_node_ids[:, el]]'

        # Compute node coordinates
        calc_node_coordinates!(node_coords, el, nodes, corners)

        # Compute Jacobian matrices
        calc_metric_terms!(jac_matrix, el, nodes, corners)

        # Compute normal vectors
        calc_normal_directions!(normal_dirs, el, nodes, corners)

        # Compute Jacobian determinants
        for i in 1:nnodes, j in 1:nnodes
            jac_det[i, j, el] = det(jac_matrix[:, :, i, j, el])
        end
    end

    # 4. Build face map (single source of truth for all faces)
    @info "Building face map..."
    face_map = build_face_map(trixi_mesh)
    determine_face_flips!(face_map, dof_id)

    # Run geometric checks if enabled
    if get(ENV, "SEAS_SEME_GEOMETRIC_CHECKS", "0") == "1"
        @info "Running geometric invariant checks..."
        try
            check_jacobian_determinant_all(jac_matrix)
            @info "  ✓ Jacobian determinants positive"

            # Check edge lengths for sample of boundary faces
            for (boundary_name, face_indices) in face_map.boundary_faces
                for idx in face_indices[1:min(5, length(face_indices))]  # Check first 5 faces
                    face_info = face_map.faces[idx]
                    face = face_info.face
                    el = face.element
                    corners = trixi_mesh.corners[:, trixi_mesh.element_node_ids[:, el]]'
                    geom = face_geometry(face, nodes, corners, jac_matrix[:, :, :, :, el])
                    check_edge_length(geom, basis.weights)
                end
            end
            @info "  ✓ Edge length checks passed"
        catch e
            @warn "Geometric check failed" exception=e
        end
    end

    # 5. Extract boundary information
    @info "Extracting boundary nodes..."

    # Material properties for impedance
    ρ = physics.density
    vs = physics.shear_velocity
    impedance_absorbing = ρ * vs
    impedance_fault = 1.0  # Fault/creep use unit impedance

    # P-wave impedance (zero for antiplane, computed for plane_strain)
    # Compute vp from poisson_ratio if plane_strain
    impedance_absorbing_p = 0.0
    if physics.formulation == :plane_strain
        μ_val = ρ * vs^2
        λ_val = 2μ_val * physics.poisson_ratio / (1 - 2 * physics.poisson_ratio)
        vp = sqrt((λ_val + 2μ_val) / ρ)
        impedance_absorbing_p = ρ * vp
    end

    # Extract fault boundary (right, upper half: 20-40 km)
    fault_ids, fault_x, fault_y, fault_mat, fault_tan, fault_nrm, fault_mat_p = get_boundary_nodes(
        trixi_mesh, node_coords, jac_matrix, basis.weights,
        impedance_fault, dof_id, :fault, nodes, face_map,
        impedance_p=0.0
    )

    # Extract creep boundary (right, lower half: 0-20 km)
    creep_ids, creep_x, creep_y, creep_mat, creep_tan, creep_nrm, creep_mat_p = get_boundary_nodes(
        trixi_mesh, node_coords, jac_matrix, basis.weights,
        impedance_fault, dof_id, :creep, nodes, face_map,
        impedance_p=0.0
    )

    # Extract absorbing boundary (left + bottom)
    absorb_ids, absorb_x, absorb_y, absorb_mat, absorb_tan, absorb_nrm, absorb_mat_p = get_boundary_nodes(
        trixi_mesh, node_coords, jac_matrix, basis.weights,
        impedance_absorbing, dof_id, :absorbing, nodes, face_map,
        impedance_p=impedance_absorbing_p
    )

    # Create BoundaryData structs
    fault_boundary = BoundaryData(
        fault_ids,
        hcat(fault_x, fault_y)',  # [2 x nnodes]
        fault_mat,
        fault_tan,
        fault_nrm,
        fault_mat_p
    )

    creep_boundary = BoundaryData(
        creep_ids,
        hcat(creep_x, creep_y)',
        creep_mat,
        creep_tan,
        creep_nrm,
        creep_mat_p
    )

    absorbing_boundary = BoundaryData(
        absorb_ids,
        hcat(absorb_x, absorb_y)',
        absorb_mat,
        absorb_tan,
        absorb_nrm,
        absorb_mat_p
    )

    boundaries = BoundaryInfo(fault_boundary, creep_boundary, absorbing_boundary)

    @info "Boundary nodes: fault=$(length(fault_ids)), creep=$(length(creep_ids)), absorbing=$(length(absorb_ids))"

    # Compute active fault mask: exclude endpoint nodes with pathological fault_mat
    # Criterion 1: nodes shared between fault and creep boundaries
    shared_creep = Set(intersect(fault_ids, creep_ids))
    # Criterion 2: nodes with fault_mat < 0.15 * median(fault_mat)
    # Catches pathological endpoint nodes (fault_mat < ~30) while keeping
    # well-conditioned interior nodes on shorter elements (fault_mat ~60-90)
    fm_sorted = sort(fault_mat)
    fm_median = fm_sorted[div(length(fm_sorted), 2)]
    fm_threshold = 0.15 * fm_median
    active_fault_mask = BitVector([
        !(fault_ids[i] in shared_creep) && fault_mat[i] >= fm_threshold
        for i in eachindex(fault_ids)
    ])
    n_excluded = count(!, active_fault_mask)
    @info "Active fault mask: $(count(active_fault_mask))/$(length(fault_ids)) active, $n_excluded excluded (threshold=$fm_threshold, median=$fm_median)"

    # Construct complete mesh
    return UnstructuredSEMesh(
        trixi_mesh,
        dof_id,
        node_coords,
        jac_matrix,
        jac_det,
        normal_dirs,
        boundaries,
        face_map,
        active_fault_mask,
        ndof,
        n_elements,
        polydeg
    )
end
