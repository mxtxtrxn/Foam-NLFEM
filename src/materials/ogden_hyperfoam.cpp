#include "foam/materials/ogden_hyperfoam.hpp"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace foam {

OgdenHyperfoam::OgdenHyperfoam(std::vector<Term> terms) : terms_(std::move(terms)) {
  for (const auto& term : terms_) {
    if (term.alpha == 0.0) {
      throw std::invalid_argument("OgdenHyperfoam: alpha_i must be nonzero");
    }
  }
}

Real OgdenHyperfoam::strain_energy(const Vector3& stretches) const {
  Real J = stretches.prod();
  if (J <= 0.0) {
    throw std::runtime_error("OgdenHyperfoam::strain_energy: J = lambda1*lambda2*lambda3 <= 0");
  }
  Real W = 0.0;
  for (const auto& term : terms_) {
    Real sum_pow = std::pow(stretches(0), term.alpha) + std::pow(stretches(1), term.alpha) +
                    std::pow(stretches(2), term.alpha) - 3.0;
    Real vol_part = (term.beta == 0.0)
                         ? -term.alpha * std::log(J)  // lim_{beta->0} (J^{-a*b}-1)/beta = -a*ln(J)
                         : (std::pow(J, -term.alpha * term.beta) - 1.0) / term.beta;
    W += (2.0 * term.mu / (term.alpha * term.alpha)) * (sum_pow + vol_part);
  }
  return W;
}

Matrix3 OgdenHyperfoam::cauchy_stress(const Matrix3& F) const {
  Real J = F.determinant();
  if (J <= 0.0) {
    throw std::runtime_error("OgdenHyperfoam::cauchy_stress: element inverted (J <= 0)");
  }

  // Spectral decomposition of the left Cauchy-Green tensor B = F F^T:
  // its eigenvectors n_k are the SPATIAL principal directions (the basis
  // in which Cauchy stress is diagonal), and its eigenvalues are
  // lambda_k^2. B (unlike F) is always symmetric positive-semidefinite,
  // so a SelfAdjointEigenSolver is both correct and numerically robust,
  // including when two or more stretches coincide.
  Eigen::SelfAdjointEigenSolver<Matrix3> eigen_solver(F * F.transpose());
  if (eigen_solver.info() != Eigen::Success) {
    throw std::runtime_error("OgdenHyperfoam::cauchy_stress: eigendecomposition of B failed");
  }
  const Vector3& lambda_squared = eigen_solver.eigenvalues();
  const Matrix3& directions = eigen_solver.eigenvectors();  // columns = n_k

  // sigma_k = sum_i (2*mu_i / (alpha_i * J)) * (lambda_k^alpha_i - J^(-alpha_i*beta_i))
  Vector3 principal_stress;
  for (int k = 0; k < 3; ++k) {
    Real lambda_k = std::sqrt(std::max(lambda_squared(k), Real(0)));
    Real sigma_k = 0.0;
    for (const auto& term : terms_) {
      sigma_k += (2.0 * term.mu / (term.alpha * J)) *
                 (std::pow(lambda_k, term.alpha) - std::pow(J, -term.alpha * term.beta));
    }
    principal_stress(k) = sigma_k;
  }

  // sigma = sum_k sigma_k * (n_k tensor n_k) = N * diag(sigma) * N^T
  return directions * principal_stress.asDiagonal() * directions.transpose();
}

Matrix6 OgdenHyperfoam::tangent_moduli(const Matrix3& F) const {
  // No analytical tangent has been derived for this model yet (see
  // docs/roadmap.md, step 1). The numerical (finite-difference) tangent
  // is a legitimate, correct-by-construction stand-in: it is exactly
  // what Newton-Raphson needs, just slower than a closed form would be.
  return numerical_tangent_moduli(*this, F);
}

}  // namespace foam
