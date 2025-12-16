"""
Discretization module: Spectral element matrices and operators.

This module provides the spatial discretization for the spectral element method:
- Mass matrices (diagonal, lumped)
- Stiffness matrices (elemental and global assembly)
- Matrix-free stiffness operator for iterative solvers
"""

# Module exports
export elemental_mass_matrix, global_mass_matrix!, build_mass_matrices
export elemental_stiffness_matrix, build_stiffness_matrices
export apply_stiffness!, StiffnessOperator, stiffness_assembly
