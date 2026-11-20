#pragma once

#include <vector>

#include "foam/materials/material_base.hpp"

namespace foam {

// Ogden hyperfoam model (the same family used by Abaqus's *HYPERFOAM):
//
//   W = sum_i (2*mu_i / alpha_i^2) * (lambda1^alpha_i + lambda2^alpha_i
//       + lambda3^alpha_i - 3 + (1/beta_i) * (J^(-alpha_i*beta_i) - 1))
//
// where lambda_k are principal stretches, J = det(F), and each term i
// contributes its own (mu_i, alpha_i, beta_i) triplet. beta_i controls
// compressibility and is what lets this model capture a near-zero
// effective Poisson's ratio typical of open-cell foams.
//
// THIS IS ROADMAP STEP 1: implement cauchy_stress()/tangent_moduli() and
// validate on make_single_hex_element() against published or analytical
// uniaxial-compression curves before touching the FEM assembly at all.
class OgdenHyperfoam : public MaterialModel {
 public:
  struct Term {
    Real mu;
    Real alpha;
    Real beta;
  };

  explicit OgdenHyperfoam(std::vector<Term> terms);

  Matrix3 cauchy_stress(const Matrix3& F) const override;
  Matrix6 tangent_moduli(const Matrix3& F) const override;

  const std::vector<Term>& terms() const { return terms_; }

 private:
  std::vector<Term> terms_;
};

}  // namespace foam
