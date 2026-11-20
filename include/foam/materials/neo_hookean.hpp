#pragma once

#include "foam/materials/material_base.hpp"

namespace foam {

// Compressible Neo-Hookean solid. Deliberately simple and well-documented
// in every textbook, so it's the material to use FIRST when testing the
// kinematics/elements/assembly/solver pipeline -- if something's broken
// in the FEM machinery, you want to find out here, not while also
// debugging the more delicate Ogden hyperfoam law.
class NeoHookean : public MaterialModel {
 public:
  NeoHookean(Real mu, Real lambda);

  Matrix3 cauchy_stress(const Matrix3& F) const override;
  Matrix6 tangent_moduli(const Matrix3& F) const override;

 private:
  Real mu_;      // shear modulus
  Real lambda_;  // Lame's first parameter
};

}  // namespace foam
