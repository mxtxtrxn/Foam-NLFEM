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

  // The LAGRANGIAN (material) tangent moduli, in Voigt notation (6x6):
  //   D_voigt = dS/dE,   S = 2nd Piola-Kirchhoff stress, E = Green
  //   -Lagrange strain, both measured in the FIXED reference
  //   configuration.
  //
  // This is deliberately NOT a spatial tangent relating a Cauchy-stress
  // rate to a rate of deformation: that version requires an objective
  // (e.g. Truesdell) stress rate, and a naive finite difference of
  // cauchy_stress() measures the PLAIN rate instead, which differs from
  // the objective one by extra convective terms (D*sigma + sigma*D -
  // tr(D)*sigma) whenever sigma != 0 -- an easy and easy-to-miss source
  // of bugs (it bit the first version of this interface, caught by the
  // cross-check test in tests/test_materials.cpp). S and E have no such
  // ambiguity: both are frame-indifferent (objective) by construction,
  // so dS/dE can be computed directly, with no rate/rotation subtlety.
  //
  // This matches a standard Total Lagrangian FEM formulation (the
  // approach planned for Quad4Element in roadmap step 2): element
  // stiffness = material part (B_L^T D B_L) + geometric part (built
  // from S directly), entirely in the reference configuration.
  virtual Matrix6 tangent_moduli(const Matrix3& F) const = 0;
};

// Generic numerical (central finite-difference) Lagrangian tangent,
// valid for ANY MaterialModel -- it only calls cauchy_stress(), nothing
// model-specific. Two uses: (a) as a material's own tangent_moduli()
// while no analytical one has been derived/verified yet (see
// OgdenHyperfoam), and (b) as an independent reference to validate an
// analytical tangent against (see the NeoHookean cross-check test).
//
// Perturbs C = F^T F directly (by +-2*dE, since E=(C-I)/2) and, for
// each perturbed C, builds a REPRESENTATIVE deformation gradient via
// the symmetric matrix square root C_perturbed^(1/2) -- i.e. a pure
// stretch, no rotation. This is valid precisely because any physically
// admissible (objective) MaterialModel must give a cauchy_stress() that
// depends on F only through C = F^T F (up to the superposed rotation
// sigma(QF) = Q sigma(F) Q^T for any rotation Q) -- so S = J F^-1 sigma
// F^-T, which strips that rotation back out, cannot depend on which F
// was used to represent a given C. Choosing the square-root (R = I) F
// is simply the most convenient representative.
//
// Voigt convention (fixed here, and binding on every later consumer of
// Matrix6, e.g. Quad4Element in roadmap step 2):
//   stress-like vector (S, or Cauchy sigma) = [s11, s22, s33, s12, s13, s23]
//   strain-like vector (E)                  = [e11, e22, e33, 2*e12, 2*e13, 2*e23]
//                                              (engineering shear strain)
// so that dot(stress_voigt, strain_voigt) equals the full tensor double
// contraction, and -- because the underlying 4th-order tensor has the
// usual minor symmetries -- Voigt matrix entries equal tensor
// components read off directly at the representative index pairs, with
// no extra bookkeeping factors.
Matrix6 numerical_tangent_moduli(const MaterialModel& material, const Matrix3& F,
                                  Real perturbation = 1e-6);

// 2nd Piola-Kirchhoff stress, derived from Cauchy stress: S = J F^-1
// sigma F^-T. Every MaterialModel only has to implement cauchy_stress(),
// but a Total Lagrangian element formulation (Quad4Element, roadmap
// step 2) works with S directly, so this conversion is exposed here
// rather than duplicated in element code.
Matrix3 second_piola_kirchhoff_stress(const MaterialModel& material, const Matrix3& F);

}  // namespace foam
