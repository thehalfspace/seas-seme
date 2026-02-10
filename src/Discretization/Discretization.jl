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

# Plane-strain exports
export elemental_stiffness_matrix_plane_strain, build_stiffness_matrices_plane_strain
export build_mass_matrices_plane_strain
export apply_stiffness_plane_strain!, StiffnessOperatorPlaneStrain
export stiffness_assembly_plane_strain
