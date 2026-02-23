"""
Discretization module: Spectral element matrices and operators.

This module provides the spatial discretization for the spectral element method:
- Mass matrices (diagonal, lumped)
- Metric weight arrays (compact per-element storage for tensor-product matvec)
- Matrix-free stiffness operator for iterative solvers
"""

# Module exports
export elemental_mass_matrix, global_mass_matrix!, build_mass_matrices
export elemental_stiffness_matrix, build_stiffness_matrices

# Metric weights (new compact storage replacing full K_el)
export MetricWeightsAntiplane, MetricWeightsPlaneStrain
export compute_metric_weights_antiplane, compute_metric_weights_plane_strain
export build_metric_weights_antiplane, build_metric_weights_plane_strain
export apply_element_stiffness_antiplane!, apply_element_stiffness_plane_strain!
export materialize_K_el_antiplane, materialize_K_el_plane_strain

# Matrix-free operators (tensor-product based)
export apply_stiffness!, StiffnessOperator, stiffness_assembly
export build_mass_matrices_plane_strain
export apply_stiffness_plane_strain!, StiffnessOperatorPlaneStrain
export stiffness_assembly_plane_strain
export elemental_stiffness_matrix_plane_strain, build_stiffness_matrices_plane_strain
