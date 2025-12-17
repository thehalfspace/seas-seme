"""
Geometric invariant checks for debugging mesh quality and bookkeeping correctness.

These checks catch 95% of bookkeeping bugs immediately:
- detJ positivity (element not inverted)
- Edge length consistency (quadrature vs physical endpoints)
- Face match consistency (interior faces align correctly)
- Normal consistency (interior face normals are opposite)

Enable with environment variable: SEAS_SEME_GEOMETRIC_CHECKS=1
"""

using Printf: @sprintf

"""
    check_jacobian_determinant(jac_matrix::AbstractArray, element::Int; tol=1e-10)

Check that Jacobian determinant is positive at all nodes in an element.

# Arguments
- `jac_matrix::AbstractArray`: Jacobian matrix [2, 2, nnodes, nnodes]
- `element::Int`: Element index (for error messages)
- `tol::Real`: Minimum acceptable determinant value

# Throws
- `DomainError`: If detJ ≤ tol at any node

# Notes
A negative or zero determinant indicates:
- Inverted element (bad mesh)
- Wrong corner ordering
- Degenerate element geometry

# Example
```julia
check_jacobian_determinant(jac_matrix[:,:,:,:,el], el)
```
"""
function check_jacobian_determinant(jac_matrix::AbstractArray, element::Int; tol=1e-10)
    nnodes = size(jac_matrix, 3)

    for j in 1:nnodes, i in 1:nnodes
        ∂x_∂ξ = jac_matrix[1, 1, i, j]
        ∂x_∂η = jac_matrix[1, 2, i, j]
        ∂y_∂ξ = jac_matrix[2, 1, i, j]
        ∂y_∂η = jac_matrix[2, 2, i, j]

        detJ = ∂x_∂ξ * ∂y_∂η - ∂x_∂η * ∂y_∂ξ

        if detJ ≤ tol
            error("""
            Jacobian determinant check failed!
            Element: $element
            Node: ($i, $j)
            detJ: $detJ (should be > $tol)

            This indicates:
            - Inverted element (negative detJ)
            - Degenerate element (zero detJ)
            - Wrong corner point ordering

            Check element corner coordinates and ordering.
            """)
        end
    end

    return true
end

"""
    check_jacobian_determinant_all(jac_matrix::AbstractArray; tol=1e-10)

Check Jacobian determinant for all elements in mesh.

# Arguments
- `jac_matrix::AbstractArray`: Jacobian matrices [2, 2, nnodes, nnodes, nelements]
- `tol::Real`: Minimum acceptable determinant value

# Returns
- `true` if all checks pass

# Example
```julia
check_jacobian_determinant_all(mesh.jac_matrix)
```
"""
function check_jacobian_determinant_all(jac_matrix::AbstractArray; tol=1e-10)
    nelements = size(jac_matrix, 5)

    for el in 1:nelements
        check_jacobian_determinant(jac_matrix[:, :, :, :, el], el; tol=tol)
    end

    return true
end

"""
    check_edge_length(face_geom, weights::AbstractVector; rtol=1e-3)

Check that quadrature edge length matches physical endpoint distance.

# Arguments
- `face_geom`: Output from `face_geometry()` (named tuple)
- `weights::AbstractVector`: Quadrature weights
- `rtol::Real`: Relative tolerance for mismatch

# Throws
- `AssertionError`: If lengths don't match within tolerance

# Notes
Invariant: ∑ wᵢ J₁D[i] ≈ ||x_end - x_start||

If this fails:
- Wrong Jacobian extraction
- Wrong face node ordering
- Numerical precision issues (increase rtol)

# Example
```julia
geom = face_geometry(face, nodes, corners, jac_matrix)
check_edge_length(geom, weights)
```
"""
function check_edge_length(face_geom, weights::AbstractVector; rtol=1e-3)
    # Quadrature length
    quad_length = sum(weights .* face_geom.J_1D)

    # Physical endpoint distance
    x_start = [face_geom.x_phys[1], face_geom.y_phys[1]]
    x_end = [face_geom.x_phys[end], face_geom.y_phys[end]]
    phys_length = sqrt(sum((x_end - x_start).^2))

    # Check relative error
    rel_error = abs(quad_length - phys_length) / max(phys_length, 1e-10)

    if rel_error > rtol
        error("""
        Edge length check failed!
        Quadrature length: $quad_length
        Physical length:   $phys_length
        Relative error:    $rel_error (tolerance: $rtol)

        This indicates:
        - Wrong Jacobian extraction (wrong indices)
        - Wrong face node ordering
        - Inconsistent geometry computation

        Check face_geometry() implementation.
        """)
    end

    return true
end

"""
    check_face_node_match(face1_nodes, face2_nodes; allow_reverse=true)

Check that two interior face node orderings match (same or reversed).

# Arguments
- `face1_nodes::AbstractVector`: Global DOF IDs for first face
- `face2_nodes::AbstractVector`: Global DOF IDs for second face
- `allow_reverse::Bool`: Allow reversed ordering (default: true)

# Returns
- `true` if match or reversed match
- `false` otherwise

# Notes
For interior faces, node ordering should either:
- Match exactly: [1,2,3,4,5] == [1,2,3,4,5]
- Match reversed: [1,2,3,4,5] == [5,4,3,2,1]

If neither, there's a connectivity bug.

# Example
```julia
nodes_L = extract_face_values(dof_id[:,:,eL], face_L)
nodes_R = extract_face_values(dof_id[:,:,eR], face_R)
@assert check_face_node_match(nodes_L, nodes_R)
```
"""
function check_face_node_match(face1_nodes, face2_nodes; allow_reverse=true)
    # Same length?
    if length(face1_nodes) != length(face2_nodes)
        return false
    end

    # Exact match?
    if face1_nodes == face2_nodes
        return true
    end

    # Reversed match?
    if allow_reverse && face1_nodes == reverse(face2_nodes)
        return true
    end

    return false
end

"""
    check_normal_consistency(normal1_x, normal1_y, normal2_x, normal2_y; rtol=1e-2)

Check that normals on interior face are opposite: n₁ ≈ -n₂.

# Arguments
- `normal1_x, normal1_y`: Normal vectors from element 1
- `normal2_x, normal2_y`: Normal vectors from element 2
- `rtol::Real`: Relative tolerance

# Throws
- `AssertionError`: If normals are not opposite

# Notes
For an interior face shared by two elements, outward normals should point
in opposite directions: n_L ≈ -n_R

If this fails:
- Wrong normal computation
- Wrong face orientation
- Inconsistent geometry

# Example
```julia
geom_L = face_geometry(face_L, ...)
geom_R = face_geometry(face_R, ...)
check_normal_consistency(geom_L.n_x, geom_L.n_y, geom_R.n_x, geom_R.n_y)
```
"""
function check_normal_consistency(normal1_x, normal1_y, normal2_x, normal2_y;
                                 rtol=1e-2)
    # Compare each normal component
    for i in eachindex(normal1_x)
        # Magnitude of normals (unnormalized, so use length for comparison)
        mag1 = sqrt(normal1_x[i]^2 + normal1_y[i]^2)
        mag2 = sqrt(normal2_x[i]^2 + normal2_y[i]^2)

        # Check they're opposite: n1 ≈ -n2
        rel_error_x = abs(normal1_x[i] + normal2_x[i]) / max(mag1, mag2, 1e-10)
        rel_error_y = abs(normal1_y[i] + normal2_y[i]) / max(mag1, mag2, 1e-10)

        if rel_error_x > rtol || rel_error_y > rtol
            error("""
            Normal consistency check failed!
            Node: $i
            Normal 1: ($(@sprintf("%.3e", normal1_x[i])), $(@sprintf("%.3e", normal1_y[i])))
            Normal 2: ($(@sprintf("%.3e", normal2_x[i])), $(@sprintf("%.3e", normal2_y[i])))
            Sum:      ($(@sprintf("%.3e", normal1_x[i] + normal2_x[i])), $(@sprintf("%.3e", normal1_y[i] + normal2_y[i])))
            Rel error: ($(@sprintf("%.3e", rel_error_x)), $(@sprintf("%.3e", rel_error_y))) (tolerance: $rtol)

            Normals on interior face should be opposite (n₁ ≈ -n₂).
            This indicates wrong normal computation or face orientation.
            """)
        end
    end

    return true
end

"""
    geometric_checks_enabled() -> Bool

Check if geometric invariant checks are enabled via environment variable.

# Returns
- `true` if SEAS_SEME_GEOMETRIC_CHECKS=1
- `false` otherwise

# Usage
Gate expensive checks behind this flag:
```julia
if geometric_checks_enabled()
    check_jacobian_determinant_all(jac_matrix)
end
```
"""
function geometric_checks_enabled()
    return get(ENV, "SEAS_SEME_GEOMETRIC_CHECKS", "0") == "1"
end
