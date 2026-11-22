#include "foam/materials/material_base.hpp"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace foam {

namespace {

// Symmetric strain tensor for Voigt index j (0..5), scaled by h. Shear
// entries are h/2 so that the *engineering* shear strain (2*epsilon_ij)
// equals h, matching the Voigt convention documented in
// material_base.hpp.
Matrix3 voigt_strain_basis(int j, Real h) {
  Matrix3 e = Matrix3::Zero();
  switch (j) {
    case 0:
      e(0, 0) = h;
      break;
    case 1:
      e(1, 1) = h;
      break;
    case 2:
      e(2, 2) = h;
      break;
    case 3:
      e(0, 1) = e(1, 0) = h / 2.0;
      break;
    case 4:
      e(0, 2) = e(2, 0) = h / 2.0;
      break;
    case 5:
      e(1, 2) = e(2, 1) = h / 2.0;
      break;
    default:
      break;
  }
  return e;
}

// Packs a symmetric tensor into Voigt order [11,22,33,12,13,23]. Used
// for both stress-like tensors (sigma, S) and, here, for C^-1-type
// tensors that appear the same way in an elasticity-tensor formula.
Eigen::Matrix<Real, 6, 1> to_voigt(const Matrix3& tensor) {
  Eigen::Matrix<Real, 6, 1> v;
  v << tensor(0, 0), tensor(1, 1), tensor(2, 2), tensor(0, 1), tensor(0, 2), tensor(1, 2);
  return v;
}

// Symmetric positive-definite matrix square root via eigendecomposition:
// sqrt(C) = V * diag(sqrt(eigenvalues)) * V^T. Used to build a
// REPRESENTATIVE (pure-stretch, zero-rotation) deformation gradient for
// a given right Cauchy-Green tensor C -- valid because 2nd
// Piola-Kirchhoff stress S depends on F only through C (by material
// objectivity), so any F with F^T F = C gives the same S.
Matrix3 symmetric_sqrt(const Matrix3& C) {
  Eigen::SelfAdjointEigenSolver<Matrix3> solver(C);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("symmetric_sqrt: eigendecomposition failed");
  }
  Vector3 sqrt_eigenvalues = solver.eigenvalues().cwiseMax(Real(0)).cwiseSqrt();
  return solver.eigenvectors() * sqrt_eigenvalues.asDiagonal() * solver.eigenvectors().transpose();
}

}  // namespace

Matrix3 second_piola_kirchhoff_stress(const MaterialModel& material, const Matrix3& F) {
  Real J = F.determinant();
  if (J <= 0.0) {
    throw std::runtime_error("second_piola_kirchhoff_stress: element inverted (J <= 0)");
  }
  Matrix3 sigma = material.cauchy_stress(F);
  Matrix3 F_inv = F.inverse();
  return J * F_inv * sigma * F_inv.transpose();
}

Matrix6 numerical_tangent_moduli(const MaterialModel& material, const Matrix3& F,
                                  Real perturbation) {
  Matrix3 C0 = F.transpose() * F;

  Matrix6 D = Matrix6::Zero();
  for (int col = 0; col < 6; ++col) {
    Matrix3 dE = voigt_strain_basis(col, perturbation);
    Matrix3 C_plus = C0 + 2.0 * dE;   // E = (C - I)/2  =>  dC = 2*dE
    Matrix3 C_minus = C0 - 2.0 * dE;

    Matrix3 F_plus = symmetric_sqrt(C_plus);
    Matrix3 F_minus = symmetric_sqrt(C_minus);

    Matrix3 S_plus = second_piola_kirchhoff_stress(material, F_plus);
    Matrix3 S_minus = second_piola_kirchhoff_stress(material, F_minus);

    Matrix3 dS = (S_plus - S_minus) / (2.0 * perturbation);
    D.col(col) = to_voigt(dS);
  }
  return D;
}

}  // namespace foam
