"""
Physics module: Material properties, friction laws, and initial conditions.

This module provides the physical models for earthquake cycle simulations:
- Rate-state friction with aging law
- Material properties (density, shear velocity, shear modulus)
- Initial stress and friction parameter distributions
"""

# Module exports
export RateStateFriction, MaterialProperties, InitialConditions
export build_initial_conditions, save_initial_conditions
export build_initial_conditions_plane_strain
