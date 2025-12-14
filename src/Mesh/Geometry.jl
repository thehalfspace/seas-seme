"""
Geometric mappings and metric terms for spectral element meshes.

This module computes coordinate transformations, Jacobian matrices, and
normal vectors for straight-sided quadrilateral elements.

Based on algorithms from Kopriva's "Implementing Spectral Methods for
Partial Differential Equations" (the "blue book").
"""

"""
    straight_side_quad_map(xi, eta, corner_points) -> (x, y)

Map from reference coordinates (ξ, η) ∈ [-1,1]² to physical coordinates (x, y).

Uses bilinear mapping for straight-sided quadrilateral elements (Alg. 95, Kopriva).

# Arguments
- `xi::Real`: Reference coordinate ξ ∈ [-1, 1]
- `eta::Real`: Reference coordinate η ∈ [-1, 1]
- `corner_points::Matrix`: Corner points [4 x 2] in physical space
  - Row order: [bottom-left, bottom-right, top-right, top-left]

# Returns
- `(x, y)`: Physical coordinates

# Notes
Corner numbering convention (counterclockwise from bottom-left):
```
4 -------- 3
|          |
|          |
1 -------- 2
```
"""
function straight_side_quad_map(xi::Real, eta::Real, corner_points::Matrix)
    x = 0.25 * (
        corner_points[1, 1] * (1.0 - xi) * (1.0 - eta) +
        corner_points[2, 1] * (1.0 + xi) * (1.0 - eta) +
        corner_points[3, 1] * (1.0 + xi) * (1.0 + eta) +
        corner_points[4, 1] * (1.0 - xi) * (1.0 + eta)
    )

    y = 0.25 * (
        corner_points[1, 2] * (1.0 - xi) * (1.0 - eta) +
        corner_points[2, 2] * (1.0 + xi) * (1.0 - eta) +
        corner_points[3, 2] * (1.0 + xi) * (1.0 + eta) +
        corner_points[4, 2] * (1.0 - xi) * (1.0 + eta)
    )

    return x, y
end

"""
    straight_side_quad_map_metrics(xi, eta, corner_points) -> (X_ξ, X_η, Y_ξ, Y_η)

Compute metric terms (Jacobian matrix components) for straight-sided mapping.

Returns the derivatives of physical coordinates with respect to reference coordinates
(Alg. 100, Kopriva).

# Arguments
- `xi::Real`: Reference coordinate ξ
- `eta::Real`: Reference coordinate η
- `corner_points::Matrix`: Corner points [4 x 2]

# Returns
- `X_ξ::Real`: ∂x/∂ξ
- `X_η::Real`: ∂x/∂η
- `Y_ξ::Real`: ∂y/∂ξ
- `Y_η::Real`: ∂y/∂η

# Notes
Jacobian matrix:
```
J = [X_ξ  X_η]
    [Y_ξ  Y_η]
```
Jacobian determinant: det(J) = X_ξ * Y_η - X_η * Y_ξ
"""
function straight_side_quad_map_metrics(xi::Real, eta::Real, corner_points::Matrix)
    X_xi = 0.25 * (
        (1.0 - eta) * (corner_points[2, 1] - corner_points[1, 1]) +
        (1.0 + eta) * (corner_points[3, 1] - corner_points[4, 1])
    )

    X_eta = 0.25 * (
        (1.0 - xi) * (corner_points[4, 1] - corner_points[1, 1]) +
        (1.0 + xi) * (corner_points[3, 1] - corner_points[2, 1])
    )

    Y_xi = 0.25 * (
        (1.0 - eta) * (corner_points[2, 2] - corner_points[1, 2]) +
        (1.0 + eta) * (corner_points[3, 2] - corner_points[4, 2])
    )

    Y_eta = 0.25 * (
        (1.0 - xi) * (corner_points[4, 2] - corner_points[1, 2]) +
        (1.0 + xi) * (corner_points[3, 2] - corner_points[2, 2])
    )

    return X_xi, X_eta, Y_xi, Y_eta
end

"""
    calc_node_coordinates!(node_coordinates, element, nodes, corners)

Compute physical (x, y) coordinates for all nodes in an element.

# Arguments
- `node_coordinates::AbstractArray{<:Any,4}`: Storage array [2, nnodes, nnodes, nelements]
  (modified in-place)
- `element::Int`: Element index
- `nodes::Vector`: Reference node locations (e.g., Gauss-Lobatto points)
- `corners::Matrix`: Element corner points [4 x 2]

# Returns
- `node_coordinates`: Modified array with physical coordinates

# Notes
- `node_coordinates[1, i, j, element]` = x-coordinate of node (i,j)
- `node_coordinates[2, i, j, element]` = y-coordinate of node (i,j)
"""
function calc_node_coordinates!(
    node_coordinates::AbstractArray{<:Any,4},
    element::Int,
    nodes::Vector,
    corners::Matrix
)
    for j in eachindex(nodes), i in eachindex(nodes)
        node_coordinates[:, i, j, element] .=
            straight_side_quad_map(nodes[i], nodes[j], corners)
    end

    return node_coordinates
end

"""
    calc_metric_terms!(jacobian_matrix, element, nodes, corners)

Compute Jacobian matrix components for all nodes in an element.

# Arguments
- `jacobian_matrix::AbstractArray`: Storage array [2, 2, nnodes, nnodes, nelements]
  (modified in-place)
- `element::Int`: Element index
- `nodes::Vector`: Reference node locations
- `corners::Matrix`: Element corner points [4 x 2]

# Returns
- `jacobian_matrix`: Modified array with metric terms

# Storage format
- `jacobian_matrix[1, 1, i, j, element]` = X_ξ
- `jacobian_matrix[1, 2, i, j, element]` = X_η
- `jacobian_matrix[2, 1, i, j, element]` = Y_ξ
- `jacobian_matrix[2, 2, i, j, element]` = Y_η
"""
function calc_metric_terms!(
    jacobian_matrix::AbstractArray,
    element::Int,
    nodes::Vector,
    corners::Matrix
)
    for j in eachindex(nodes), i in eachindex(nodes)
        (
            jacobian_matrix[1, 1, i, j, element],
            jacobian_matrix[1, 2, i, j, element],
            jacobian_matrix[2, 1, i, j, element],
            jacobian_matrix[2, 2, i, j, element]
        ) = straight_side_quad_map_metrics(nodes[i], nodes[j], corners)
    end

    return jacobian_matrix
end

"""
    calc_normal_directions!(normal_directions, element, nodes, corners)

Compute outward-pointing normal vectors on element boundaries.

# Arguments
- `normal_directions::AbstractArray`: Storage array [2, nnodes, 4, nelements]
  (modified in-place)
- `element::Int`: Element index
- `nodes::Vector`: Reference node locations
- `corners::Matrix`: Element corner points [4 x 2]

# Returns
- `normal_directions`: Modified array with normal vectors

# Notes
- Surface numbering (counterclockwise from bottom):
  - Side 1: Bottom (η = -1)
  - Side 2: Right (ξ = 1)
  - Side 3: Top (η = 1)
  - Side 4: Left (ξ = -1)
- Normal vectors are NOT normalized (normalization done during flux computation)
- `normal_directions[1, i, side, element]` = x-component of normal
- `normal_directions[2, i, side, element]` = y-component of normal
"""
function calc_normal_directions!(
    normal_directions::AbstractArray,
    element::Int,
    nodes::Vector,
    corners::Matrix
)
    # Sides 2 (right) and 4 (left)
    for j in eachindex(nodes)
        # Side 2: ξ = 1 (right edge)
        X_xi, X_eta, Y_xi, Y_eta = straight_side_quad_map_metrics(1.0, nodes[j], corners)
        Jtemp = X_xi * Y_eta - X_eta * Y_xi
        normal_directions[1, j, 2, element] = sign(Jtemp) * Y_eta
        normal_directions[2, j, 2, element] = sign(Jtemp) * (-X_eta)

        # Side 4: ξ = -1 (left edge)
        X_xi, X_eta, Y_xi, Y_eta = straight_side_quad_map_metrics(-1.0, nodes[j], corners)
        Jtemp = X_xi * Y_eta - X_eta * Y_xi
        normal_directions[1, j, 4, element] = -sign(Jtemp) * Y_eta
        normal_directions[2, j, 4, element] = -sign(Jtemp) * (-X_eta)
    end

    # Sides 1 (bottom) and 3 (top)
    for i in eachindex(nodes)
        # Side 1: η = -1 (bottom edge)
        X_xi, X_eta, Y_xi, Y_eta = straight_side_quad_map_metrics(nodes[i], -1.0, corners)
        Jtemp = X_xi * Y_eta - X_eta * Y_xi
        normal_directions[1, i, 1, element] = -sign(Jtemp) * (-Y_xi)
        normal_directions[2, i, 1, element] = -sign(Jtemp) * X_xi

        # Side 3: η = 1 (top edge)
        X_xi, X_eta, Y_xi, Y_eta = straight_side_quad_map_metrics(nodes[i], 1.0, corners)
        Jtemp = X_xi * Y_eta - X_eta * Y_xi
        normal_directions[1, i, 3, element] = sign(Jtemp) * (-Y_xi)
        normal_directions[2, i, 3, element] = sign(Jtemp) * X_xi
    end

    return normal_directions
end
