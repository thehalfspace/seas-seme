"""
Matrix-free stiffness operator for iterative solvers.

Implements efficient element-by-element matrix-vector products without
assembling the global stiffness matrix.
"""

"""
    apply_stiffness!(y, x, K_el, dof_id, n_elements)

Apply stiffness operator y = K*x using element-by-element operations (matrix-free).

# Arguments
- `y::Vector`: Output vector (modified in-place, zeroed first)
- `x::Vector`: Input vector
- `K_el::Array{T,3}`: Elemental stiffness matrices [nnodes², nnodes², n_elements]
- `dof_id::Array{Int,3}`: Connectivity matrix [nnodes, nnodes, n_elements]
- `n_elements::Int`: Number of elements

# Returns
- `y::Vector`: Result of K*x

# Algorithm
```
y .= 0
for each element:
    extract local x values
    compute y_local = K_el * x_local
    accumulate y_local into global y
```

# Notes
- Avoids assembling global sparse matrix (major memory savings for high-order p)
- Element-local operations have good cache locality
- Can be parallelized with mesh coloring
"""
function apply_stiffness!(
    y::Vector{T},
    x::Vector{T},
    K_el::Array{T,3},
    dof_id::Array{Int,3},
    n_elements::Int
) where {T<:AbstractFloat}
    fill!(y, zero(T))  # Zero output vector

    for el in 1:n_elements
        # Get local DOF indices for this element
        local_idx = dof_id[:, :, el]

        # Extract local vector (gather)
        x_local = x[local_idx]

        # Local matrix-vector product
        y_local = K_el[:, :, el] * x_local[:]

        # Accumulate into global vector (scatter-add)
        y[local_idx] .+= reshape(y_local, size(local_idx))
    end

    return y
end

"""
    StiffnessOperator{T} <: AbstractMatrix{T}

Matrix-free stiffness operator for use with iterative linear solvers.

# Fields
- `K_el::Array{T,3}`: Elemental stiffness matrices [nnodes², nnodes², n_elements]
- `dof_id::Array{Int,3}`: Connectivity matrix [nnodes, nnodes, n_elements]
- `n_elements::Int`: Number of elements
- `ndof::Int`: Total number of global DOFs
- `fltni::Vector{Int}`: Non-fault/creep DOF indices (interior DOFs)
- `x_full::Vector{T}`: Workspace for full DOF vector (pre-allocated)
- `y_full::Vector{T}`: Workspace for full DOF vector (pre-allocated)

# Usage with Iterative Solvers
This operator can be used with any iterative solver from IterativeSolvers.jl
or Krylov.jl that accepts `AbstractMatrix`:

```julia
using IterativeSolvers

# Create operator
K_op = StiffnessOperator(K_el, dof_id, n_elements, ndof, fltni)

# Solve K*u = -rhs with CG
u = cg(K_op, -rhs, Pl=preconditioner)
```

# Notes
- Implements `LinearAlgebra.mul!(y, A, x)` for matrix-vector products
- Only operates on non-fault DOFs (fault DOFs are prescribed)
- Pre-allocated workspaces eliminate allocations in iterative loop
- Compatible with AMG preconditioners
"""
struct StiffnessOperator{T<:AbstractFloat} <: AbstractMatrix{T}
    K_el::Array{T,3}
    dof_id::Array{Int,3}
    n_elements::Int
    ndof::Int
    fltni::Vector{Int}
    x_full::Vector{T}
    y_full::Vector{T}
end

"""
    StiffnessOperator(K_el, dof_id, n_elements, ndof, fltni)

Construct matrix-free stiffness operator with pre-allocated workspaces.
"""
function StiffnessOperator(
    K_el::Array{T,3},
    dof_id::Array{Int,3},
    n_elements::Int,
    ndof::Int,
    fltni::Vector{Int}
) where {T<:AbstractFloat}
    # Pre-allocate workspace arrays
    x_full = zeros(T, ndof)
    y_full = zeros(T, ndof)

    return StiffnessOperator(K_el, dof_id, n_elements, ndof, fltni, x_full, y_full)
end

"""
    LinearAlgebra.mul!(y, A::StiffnessOperator, x)

Matrix-vector product y = A*x for matrix-free stiffness operator.

# Arguments
- `y::Vector`: Output vector (non-fault DOFs only)
- `A::StiffnessOperator`: Matrix-free operator
- `x::Vector`: Input vector (non-fault DOFs only)

# Returns
- `y::Vector`: Result (modified in-place)

# Algorithm
1. Expand x to full DOF vector (zero on fault/creep boundaries)
2. Apply stiffness element-by-element
3. Extract non-fault DOFs from result

# Notes
- Uses pre-allocated workspaces to avoid allocations
- All operations are in-place for performance
"""
function LinearAlgebra.mul!(y::Vector, A::StiffnessOperator{T}, x::Vector) where {T}
    # Expand x to full DOF space (zero on fault/creep)
    fill!(A.x_full, zero(T))
    A.x_full[A.fltni] .= x

    # Apply stiffness element-by-element (matrix-free!)
    apply_stiffness!(A.y_full, A.x_full, A.K_el, A.dof_id, A.n_elements)

    # Extract non-fault DOFs
    y .= A.y_full[A.fltni]

    return y
end

# Implement required LinearAlgebra interface
Base.size(A::StiffnessOperator) = (length(A.fltni), length(A.fltni))
Base.size(A::StiffnessOperator, d::Int) = size(A)[d]
Base.eltype(A::StiffnessOperator{T}) where {T} = T

"""
    stiffness_assembly(K_el, dof_id) -> SparseMatrixCSC

Assemble global sparse stiffness matrix from elemental matrices.

# Arguments
- `K_el::Array{T,3}`: Elemental stiffness matrices
- `dof_id::Array{Int,3}`: Connectivity matrix

# Returns
- `SparseMatrixCSC`: Assembled global stiffness matrix

# Notes
- Used ONLY for building AMG preconditioner (one-time cost)
- Main solve uses matrix-free operator for efficiency
- Assembles by collecting (I, J, V) triplets then constructing sparse matrix
"""
function stiffness_assembly(K_el::Array{T,3}, dof_id::Array{Int,3}) where {T}
    # Preallocate triplet vectors
    I = Vector{Int}(undef, length(K_el))
    J = Vector{Int}(undef, length(K_el))
    V = Vector{T}(undef, length(K_el))

    Nel = size(dof_id, 3)  # Number of elements
    ndof = maximum(dof_id)

    ct = 1
    for el in 1:Nel
        v = view(dof_id, :, :, el)
        for j in 1:length(v), i in 1:length(v)
            I[ct] = v[i]
            J[ct] = v[j]
            V[ct] = K_el[i, j, el]
            ct += 1
        end
    end

    # Construct sparse matrix (sum duplicates with '+' combiner)
    return sparse(I, J, V, ndof, ndof, +)
end
