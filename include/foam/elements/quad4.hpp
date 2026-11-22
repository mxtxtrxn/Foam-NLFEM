#pragma once

#include <array>

#include "foam/elements/element_base.hpp"

namespace foam {

// 4-node bilinear quadrilateral, plane-strain, Total Lagrangian
// formulation (consistent with MaterialModel::tangent_moduli() being
// the Lagrangian dS/dE -- see docs/architecture.md). 2x2 Gauss
// quadrature. A good second step after the single Hex8 element: small
// enough to reason about by hand, general enough to build a real mesh
// with (roadmap step 2: structured grid, compression specimen).
//
// Node ordering: bottom-left, bottom-right, top-right, top-left (CCW),
// matching the standard bilinear shape functions N1..N4 at natural
// coordinates (-1,-1), (1,-1), (1,1), (-1,1). make_structured_grid_2d()
// generates meshes with this same ordering.
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

 private:
  // Everything needed at one Gauss point to build either the internal
  // force or the tangent stiffness contribution there. Computed once
  // per Gauss point and shared between internal_force() and
  // tangent_stiffness() so the Jacobian/shape-gradient/F/B_L algebra --
  // the easiest part of this element to get subtly wrong -- exists in
  // exactly one place.
  struct GaussPointData {
    Real detJ0;                        // reference-configuration Jacobian determinant
    Eigen::Matrix<Real, 4, 2> dN_dX;    // reference shape-function gradients, one row per node
    Matrix3 F;                         // deformation gradient (plane-strain embedded)
    Eigen::Matrix<Real, 3, 8> B_L;      // nonlinear strain-displacement matrix, rows [E11,E22,2E12]
  };

  std::array<GaussPointData, 4> compute_gauss_point_data(const MatrixX& reference_positions,
                                                          const MatrixX& nodal_displacements) const;
};

}  // namespace foam
