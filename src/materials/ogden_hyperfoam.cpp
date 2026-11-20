#include "foam/materials/ogden_hyperfoam.hpp"

#include <stdexcept>

namespace foam {

OgdenHyperfoam::OgdenHyperfoam(std::vector<Term> terms) : terms_(std::move(terms)) {}

Matrix3 OgdenHyperfoam::cauchy_stress(const Matrix3& /*F*/) const {
  // TODO(roadmap step 1):
  //   1. Eigen-decompose F (or C = F^T F) to get principal stretches
  //      lambda1, lambda2, lambda3 and principal directions.
  //   2. Evaluate dW/d(lambda_k) per term, sum over terms_.
  //   3. Convert principal Kirchhoff/Cauchy stress back to the
  //      Cartesian basis using the principal directions.
  throw std::logic_error("OgdenHyperfoam::cauchy_stress: not implemented yet");
}

Matrix6 OgdenHyperfoam::tangent_moduli(const Matrix3& /*F*/) const {
  // TODO(roadmap step 1): consistent tangent. Start by validating
  // cauchy_stress() against finite differences before deriving this
  // analytically -- a numerically-differentiated tangent is a
  // legitimate (if slower) stand-in for early single-element tests.
  throw std::logic_error("OgdenHyperfoam::tangent_moduli: not implemented yet");
}

}  // namespace foam
