#pragma once

#include "foam/core/types.hpp"

namespace foam {

// Interface every constitutive model must satisfy. Everything above this
// layer (elements, assembly, solvers) talks to materials ONLY through
// this interface -- so a new model (visco-hyperelastic, micromechanical
// RVE surrogate, neural-network surrogate, ...) can be dropped in without
// touching the FEM machinery.
class MaterialModel {
 public:
  virtual ~MaterialModel() = default;

  // Cauchy stress (current configuration) for a given deformation
  // gradient F.
  virtual Matrix3 cauchy_stress(const Matrix3& F) const = 0;

  // Spatial tangent moduli in Voigt notation (6x6), i.e. d(sigma)/d(epsilon)
  // consistently linearized for use in Newton-Raphson. This is the
  // expensive, error-prone part to derive by hand for an Ogden-type
  // energy -- a good candidate to cross-check with finite differences
  // or symbolic/automatic differentiation in a unit test.
  virtual Matrix6 tangent_moduli(const Matrix3& F) const = 0;
};

}  // namespace foam
