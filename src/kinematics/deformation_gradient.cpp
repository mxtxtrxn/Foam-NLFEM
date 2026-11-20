#include "foam/kinematics/deformation_gradient.hpp"

#include <stdexcept>

namespace foam {

Matrix3 compute_deformation_gradient(const MatrixX& /*shape_gradients*/,
                                      const MatrixX& /*nodal_displacements*/) {
  // TODO(roadmap step 1): F = I + sum_a (u_a outer grad(N_a)).
  throw std::logic_error("compute_deformation_gradient: not implemented yet");
}

Real jacobian_determinant(const Matrix3& F) { return F.determinant(); }

}  // namespace foam
