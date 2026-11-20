#pragma once

#include "foam/core/types.hpp"
#include "foam/materials/material_base.hpp"

namespace foam {

// A finite-element FORMULATION: shape functions, quadrature, and the
// element-level residual (internal force) and tangent stiffness. This is
// distinct from ElementConnectivity (foam/mesh), which only records mesh
// topology. One FiniteElement instance is typically shared across all
// elements of the same type/order; it is stateless and takes nodal
// positions/displacements/a MaterialModel as arguments.
class FiniteElement {
 public:
  virtual ~FiniteElement() = default;

  virtual int num_nodes() const = 0;
  virtual int num_dofs_per_node() const = 0;

  // Element-level internal force vector (residual contribution),
  // computed by numerically integrating the material's Cauchy stress
  // over the element.
  virtual VectorX internal_force(const MatrixX& reference_positions,
                                  const MatrixX& nodal_displacements,
                                  const MaterialModel& material) const = 0;

  // Element-level tangent stiffness (material + geometric parts),
  // consistent with internal_force for quadratic Newton convergence.
  virtual MatrixX tangent_stiffness(const MatrixX& reference_positions,
                                     const MatrixX& nodal_displacements,
                                     const MaterialModel& material) const = 0;
};

}  // namespace foam
