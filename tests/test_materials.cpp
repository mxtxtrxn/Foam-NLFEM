#include <cmath>

#include "doctest.h"
#include "foam/materials/neo_hookean.hpp"
#include "foam/materials/ogden_hyperfoam.hpp"

// ------------------------------- NeoHookean -------------------------------

TEST_CASE("NeoHookean: undeformed state (F = I) gives zero stress") {
  foam::NeoHookean material(/*mu=*/1.0, /*lambda=*/1.0);
  foam::Matrix3 F = foam::Matrix3::Identity();
  foam::Matrix3 sigma = material.cauchy_stress(F);
  CHECK(sigma.norm() == doctest::Approx(0.0).epsilon(1e-12));
}

TEST_CASE("NeoHookean: rejects an inverted element (J <= 0)") {
  foam::NeoHookean material(1.0, 1.0);
  foam::Matrix3 F = foam::Matrix3::Zero();  // det(F) = 0
  CHECK_THROWS_AS(material.cauchy_stress(F), std::runtime_error);
}

TEST_CASE("NeoHookean: closed-form tangent matches the generic numerical tangent") {
  foam::NeoHookean material(1.5e5, 2.0e5);  // Pa; arbitrary but realistic order of magnitude
  foam::Matrix3 F;
  F << 1.10, 0.05, 0.00, 0.02, 0.95, 0.00, 0.00, 0.00, 1.00;

  foam::Matrix6 C_analytical = material.tangent_moduli(F);
  foam::Matrix6 C_numerical = foam::numerical_tangent_moduli(material, F, 1e-6);

  double scale = C_analytical.cwiseAbs().maxCoeff();
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      CHECK(std::abs(C_analytical(i, j) - C_numerical(i, j)) < scale * 1e-4);
    }
  }
}

// ----------------------------- OgdenHyperfoam -----------------------------

TEST_CASE("OgdenHyperfoam: undeformed state (F = I) gives zero stress") {
  foam::OgdenHyperfoam material({foam::OgdenHyperfoam::Term{80000.0, 13.0, 0.4}});
  foam::Matrix3 F = foam::Matrix3::Identity();
  CHECK(material.cauchy_stress(F).norm() == doctest::Approx(0.0).epsilon(1e-9));
}

TEST_CASE("OgdenHyperfoam: rejects an inverted element (J <= 0)") {
  foam::OgdenHyperfoam material({foam::OgdenHyperfoam::Term{80000.0, 13.0, 0.4}});
  foam::Matrix3 F = foam::Matrix3::Zero();
  CHECK_THROWS_AS(material.cauchy_stress(F), std::runtime_error);
}

TEST_CASE(
    "OgdenHyperfoam: closed-form principal stress matches finite-difference of the strain "
    "energy (independent check of the formula, not just the implementation)") {
  foam::OgdenHyperfoam material({foam::OgdenHyperfoam::Term{80000.0, 13.0, 0.4}});

  foam::Vector3 stretches(0.8, 1.05, 1.15);  // distinct values -> no eigenvector ambiguity
  foam::Real J = stretches.prod();

  double h = 1e-6;
  foam::Vector3 expected_sigma;
  for (int k = 0; k < 3; ++k) {
    foam::Vector3 plus = stretches;
    foam::Vector3 minus = stretches;
    plus(k) += h;
    minus(k) -= h;
    double dW_dlambda = (material.strain_energy(plus) - material.strain_energy(minus)) / (2.0 * h);
    expected_sigma(k) = (stretches(k) / J) * dW_dlambda;
  }

  foam::Matrix3 F = foam::Matrix3::Zero();
  F(0, 0) = stretches(0);
  F(1, 1) = stretches(1);
  F(2, 2) = stretches(2);
  foam::Matrix3 sigma = material.cauchy_stress(F);

  CHECK(sigma(0, 0) == doctest::Approx(expected_sigma(0)).epsilon(1e-4));
  CHECK(sigma(1, 1) == doctest::Approx(expected_sigma(1)).epsilon(1e-4));
  CHECK(sigma(2, 2) == doctest::Approx(expected_sigma(2)).epsilon(1e-4));
  CHECK(std::abs(sigma(0, 1)) < 1.0);
  CHECK(std::abs(sigma(0, 2)) < 1.0);
  CHECK(std::abs(sigma(1, 2)) < 1.0);
}

TEST_CASE("OgdenHyperfoam: objective under a superposed rigid rotation (sigma' = R sigma R^T)") {
  foam::OgdenHyperfoam material({foam::OgdenHyperfoam::Term{80000.0, 13.0, 0.4}});

  foam::Matrix3 F;
  F << 0.85, 0.05, 0.00, 0.02, 1.05, 0.03, 0.00, 0.01, 1.10;

  double theta = 0.4;
  foam::Matrix3 Rz;
  Rz << std::cos(theta), -std::sin(theta), 0, std::sin(theta), std::cos(theta), 0, 0, 0, 1;
  double phi = 0.3;
  foam::Matrix3 Rx;
  Rx << 1, 0, 0, 0, std::cos(phi), -std::sin(phi), 0, std::sin(phi), std::cos(phi);
  foam::Matrix3 R = Rx * Rz;  // a general, non-axis-aligned rotation

  foam::Matrix3 sigma = material.cauchy_stress(F);
  foam::Matrix3 sigma_rotated = material.cauchy_stress(R * F);
  foam::Matrix3 expected = R * sigma * R.transpose();

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      CHECK(std::abs(sigma_rotated(i, j) - expected(i, j)) < 1e-6);
    }
  }
}

TEST_CASE("OgdenHyperfoam: beta=0 gives exactly zero lateral stress with no lateral stretch") {
  // This is the well-known "zero effective Poisson's ratio" property of
  // beta=0 terms: an unconfined compression specimen needs NO lateral
  // bulging/contraction to satisfy the free-surface (zero lateral
  // stress) condition.
  foam::OgdenHyperfoam material({foam::OgdenHyperfoam::Term{80000.0, 13.0, 0.0}});
  foam::Matrix3 F = foam::Matrix3::Zero();
  F(0, 0) = 0.6;  // substantial axial compression
  F(1, 1) = 1.0;
  F(2, 2) = 1.0;

  foam::Matrix3 sigma = material.cauchy_stress(F);
  CHECK(std::abs(sigma(1, 1)) < 1e-6);
  CHECK(std::abs(sigma(2, 2)) < 1e-6);
  CHECK(sigma(0, 0) < 0.0);  // compressive
}

TEST_CASE("second_piola_kirchhoff_stress: round-trips back to cauchy_stress via sigma = "
          "(1/J) F S F^T") {
  foam::OgdenHyperfoam material({foam::OgdenHyperfoam::Term{80000.0, 13.0, 0.4}});
  foam::Matrix3 F;
  F << 0.85, 0.05, 0.00, 0.02, 1.05, 0.03, 0.00, 0.01, 1.10;

  foam::Matrix3 sigma = material.cauchy_stress(F);
  foam::Matrix3 S = foam::second_piola_kirchhoff_stress(material, F);
  foam::Real J = F.determinant();
  foam::Matrix3 sigma_roundtrip = (1.0 / J) * F * S * F.transpose();

  CHECK((sigma - sigma_roundtrip).norm() < 1e-8);
}

TEST_CASE("OgdenHyperfoam: numerical tangent modulus is symmetric (hyperelastic consistency)") {
  foam::OgdenHyperfoam material({foam::OgdenHyperfoam::Term{80000.0, 13.0, 0.3}});
  foam::Matrix3 F = foam::Matrix3::Zero();
  F(0, 0) = 0.85;
  F(1, 1) = 1.02;
  F(2, 2) = 1.02;

  foam::Matrix6 C = material.tangent_moduli(F);
  double scale = C.cwiseAbs().maxCoeff();
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      CHECK(std::abs(C(i, j) - C(j, i)) < scale * 1e-4);
    }
  }
}
