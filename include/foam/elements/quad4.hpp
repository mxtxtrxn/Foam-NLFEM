#pragma once

#include "foam/elements/element_base.hpp"

namespace foam {

// 4-node bilinear quadrilateral, plane-strain by default. A good second
// step after the single Hex8 element: small enough to reason about by
// hand (2x2 Gauss quadrature), general enough to build a real mesh with
// (roadmap step 2: structured grid, compression specimen).
class Quad4Element : public FiniteElement {
 public:
  int num_nodes() const override { return 4; }
  int num_dofs_per_node() const override { return 2; }

  VectorX internal_force(const MatrixX& reference_positions,
                          const MatrixX& nodal_displacements,
                          const MaterialModel& material) const override;

  MatrixX tangent_stiffness(const MatrixX& reference_positions,
                             const MatrixX& nodal_displacements,
                             const MaterialModel& material) const override;
};

}  // namespace foam
