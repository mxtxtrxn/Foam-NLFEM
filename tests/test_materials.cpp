#include "doctest.h"
#include "foam/materials/neo_hookean.hpp"
#include "foam/materials/ogden_hyperfoam.hpp"

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

// --- Stub-tracking tests -------------------------------------------------
// These document CURRENT behavior (roadmap step 1 not done yet). Once
// OgdenHyperfoam::cauchy_stress/tangent_moduli are implemented for real,
// delete these two tests and replace them with ones that check actual
// numbers against a published/analytical uniaxial-compression curve.

TEST_CASE("OgdenHyperfoam::cauchy_stress is a documented not-yet-implemented stub") {
  foam::OgdenHyperfoam material({foam::OgdenHyperfoam::Term{1.0, 2.0, 0.3}});
  foam::Matrix3 F = foam::Matrix3::Identity();
  CHECK_THROWS_AS(material.cauchy_stress(F), std::logic_error);
}

TEST_CASE("OgdenHyperfoam::tangent_moduli is a documented not-yet-implemented stub") {
  foam::OgdenHyperfoam material({foam::OgdenHyperfoam::Term{1.0, 2.0, 0.3}});
  foam::Matrix3 F = foam::Matrix3::Identity();
  CHECK_THROWS_AS(material.tangent_moduli(F), std::logic_error);
}
