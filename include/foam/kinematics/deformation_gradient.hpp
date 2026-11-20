#pragma once

#include "foam/core/types.hpp"

namespace foam {

// F = I + grad(u). Computed from an element's shape-function gradients
// (in the reference configuration) and the nodal displacements.
// `shape_gradients` has one row per node, giving d(N_a)/dX in the
// reference configuration; `nodal_displacements` has one row per node.
//
// TODO(roadmap step 1): implement for the Hex8 single-element case first;
// generalize to other element types afterwards.
Matrix3 compute_deformation_gradient(const MatrixX& shape_gradients,
                                      const MatrixX& nodal_displacements);

// det(F): local volume ratio. Foams are highly compressible, so this is
// checked often (e.g. to flag element inversion, det(F) <= 0).
Real jacobian_determinant(const Matrix3& F);

}  // namespace foam
