#include "foam/kinematics/strain_measures.hpp"

namespace foam {

Invariants principal_invariants(const Matrix3& C) {
  Real trC = C.trace();
  Real trC2 = (C * C).trace();
  Invariants inv;
  inv.I1 = trC;
  inv.I2 = 0.5 * (trC * trC - trC2);
  inv.I3 = C.determinant();
  return inv;
}

}  // namespace foam
