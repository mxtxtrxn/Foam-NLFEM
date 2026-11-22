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

Matrix6 NeoHookean::tangent_moduli(const Matrix3& F) const {
  // Lagrangian (material) tangent D = dS/dE, derived directly from
  // S = mu*(I - C^-1) + lambda*ln(J)*C^-1 via dS/dC (using the standard
  // identities d(C^-1)_IJ/dC_KL = -1/2(Cinv_IK Cinv_JL + Cinv_IL Cinv_JK)
  // and d(ln J)/dC = (1/2) C^-1), then D = 2*dS/dC:
  //
  //   D_IJKL = lambda*Cinv_IJ*Cinv_KL
  //          + (mu - lambda*ln(J)) * (Cinv_IK*Cinv_JL + Cinv_IL*Cinv_JK)
  //
  // Sanity check (worth keeping in mind when reading this): at C = I
  // (undeformed), Cinv = I and ln(J) = 0, so this reduces to
  // D_IJKL = lambda*d_IJ*d_KL + mu*(d_IK*d_JL + d_IL*d_JK), i.e. exactly
  // the standard linear-isotropic-elasticity tensor -- as it must, since
  // any hyperelastic model should linearize to ordinary linear
  // elasticity at zero strain.
  Real J = F.determinant();
  if (J <= 0.0) {
    throw std::runtime_error("NeoHookean::tangent_moduli: element inverted (J <= 0)");
  }
  Matrix3 C = F.transpose() * F;
  Matrix3 Cinv = C.inverse();
  Real coef2 = mu_ - lambda_ * std::log(J);

  // Representative (i,j) tensor index pair for each Voigt index 0..5.
  static constexpr int idx_i[6] = {0, 1, 2, 0, 0, 1};
  static constexpr int idx_j[6] = {0, 1, 2, 1, 2, 2};

  Matrix6 D;
  for (int A = 0; A < 6; ++A) {
    int i = idx_i[A];
    int j = idx_j[A];
    for (int B = 0; B < 6; ++B) {
      int k = idx_i[B];
      int l = idx_j[B];
      D(A, B) = lambda_ * Cinv(i, j) * Cinv(k, l) +
                coef2 * (Cinv(i, k) * Cinv(j, l) + Cinv(i, l) * Cinv(j, k));
    }
  }
  return D;
}

}  // namespace foam
