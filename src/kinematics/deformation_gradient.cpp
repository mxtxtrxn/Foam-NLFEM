#include "foam/kinematics/deformation_gradient.hpp"

#include <stdexcept>

namespace foam {

Matrix3 compute_deformation_gradient(const MatrixX& shape_gradients,
                                      const MatrixX& nodal_displacements) {
  if (shape_gradients.rows() != nodal_displacements.rows() ||
      shape_gradients.cols() != nodal_displacements.cols()) {
    throw std::invalid_argument(
        "compute_deformation_gradient: shape_gradients and nodal_displacements must have "
        "matching dimensions (one row per node, one column per spatial dimension)");
  }
  int num_nodes = static_cast<int>(shape_gradients.rows());
  int dim = static_cast<int>(shape_gradients.cols());
  if (dim != 2 && dim != 3) {
    throw std::invalid_argument("compute_deformation_gradient: dim must be 2 or 3");
  }

  // F = I + sum_a (u_a outer grad(N_a)). For dim=2, only F's top-left
  // 2x2 block is touched -- F(2,2) stays at its identity value of 1,
  // and the rest of row/column 2 stays 0. That is exactly the standard
  // PLANE STRAIN embedding (no out-of-plane displacement, no
  // z-dependence), which is what Quad4Element relies on.
  Matrix3 F = identity3();
  for (int a = 0; a < num_nodes; ++a) {
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < dim; ++j) {
        F(i, j) += nodal_displacements(a, i) * shape_gradients(a, j);
      }
    }
  }
  return F;
}

Real jacobian_determinant(const Matrix3& F) { return F.determinant(); }

}  // namespace foam
