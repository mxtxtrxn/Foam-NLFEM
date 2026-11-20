#include "foam/materials/neo_hookean.hpp"

#include <cmath>
#include <stdexcept>

namespace foam {

NeoHookean::NeoHookean(Real mu, Real lambda) : mu_(mu), lambda_(lambda) {}

Matrix3 NeoHookean::cauchy_stress(const Matrix3& F) const {
  // Standard compressible Neo-Hookean:
  //   sigma = (mu/J) (B - I) + (lambda * ln(J) / J) I,   B = F F^T, J = det(F)
  Real J = F.determinant();
  if (J <= 0.0) {
    throw std::runtime_error("NeoHookean::cauchy_stress: element inverted (J <= 0)");
  }
  Matrix3 B = F * F.transpose();
  return (mu_ / J) * (B - identity3()) + (lambda_ * std::log(J) / J) * identity3();
}

Matrix6 NeoHookean::tangent_moduli(const Matrix3& /*F*/) const {
  // TODO(roadmap step 1, "known-good" material): the spatial tangent for
  // this model has a standard closed form,
  //   c_ijkl = (lambda/J) d_ij d_kl + [(mu - lambda*ln(J))/J] (d_ik d_jl + d_il d_jk)
  // Implement it in Voigt notation and cross-check it against a finite
  // -difference derivative of cauchy_stress() in a unit test before
  // trusting it inside Newton-Raphson.
  throw std::logic_error("NeoHookean::tangent_moduli: not implemented yet");
}

}  // namespace foam
