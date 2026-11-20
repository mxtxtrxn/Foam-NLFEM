#pragma once

#include "foam/core/types.hpp"

namespace foam {

// Right Cauchy-Green tensor: C = F^T F.
inline Matrix3 right_cauchy_green(const Matrix3& F) { return F.transpose() * F; }

// Green-Lagrange strain: E = 1/2 (C - I).
inline Matrix3 green_lagrange_strain(const Matrix3& F) {
  return 0.5 * (right_cauchy_green(F) - identity3());
}

// Principal invariants of C (I1, I2, I3), used directly by Ogden-type
// strain-energy functions expressed in terms of principal stretches,
// and by many isotropic hyperelastic models in general.
struct Invariants {
  Real I1;
  Real I2;
  Real I3;
};

Invariants principal_invariants(const Matrix3& C);

}  // namespace foam
