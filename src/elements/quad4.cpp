#include "foam/elements/quad4.hpp"

#include <stdexcept>

namespace foam {

VectorX Quad4Element::internal_force(const MatrixX& /*reference_positions*/,
                                      const MatrixX& /*nodal_displacements*/,
                                      const MaterialModel& /*material*/) const {
  // TODO(roadmap step 2):
  //   for each of the 4 Gauss points:
  //     - evaluate bilinear shape function gradients at that point
  //     - compute F there via compute_deformation_gradient()
  //     - get sigma = material.cauchy_stress(F)
  //     - accumulate B^T * sigma * detJ * weight into the residual
  throw std::logic_error("Quad4Element::internal_force: not implemented yet");
}

MatrixX Quad4Element::tangent_stiffness(const MatrixX& /*reference_positions*/,
                                         const MatrixX& /*nodal_displacements*/,
                                         const MaterialModel& /*material*/) const {
  // TODO(roadmap step 2): material tangent (via material.tangent_moduli())
  // plus the geometric stiffness term from the stress already present.
  throw std::logic_error("Quad4Element::tangent_stiffness: not implemented yet");
}

}  // namespace foam
