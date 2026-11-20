#include "doctest.h"
#include "foam/kinematics/deformation_gradient.hpp"
#include "foam/kinematics/strain_measures.hpp"

TEST_CASE("principal_invariants of C = I (undeformed state)") {
  foam::Matrix3 C = foam::Matrix3::Identity();
  auto inv = foam::principal_invariants(C);
  CHECK(inv.I1 == doctest::Approx(3.0));
  CHECK(inv.I2 == doctest::Approx(3.0));
  CHECK(inv.I3 == doctest::Approx(1.0));
}

TEST_CASE("principal_invariants of a uniform dilation C = s^2 I") {
  double s = 2.0;
  foam::Matrix3 C = (s * s) * foam::Matrix3::Identity();
  auto inv = foam::principal_invariants(C);
  CHECK(inv.I1 == doctest::Approx(3.0 * s * s));
  CHECK(inv.I3 == doctest::Approx(s * s * s * s * s * s));  // (s^2)^3
}

TEST_CASE("jacobian_determinant matches Eigen's determinant") {
  foam::Matrix3 F = 2.0 * foam::Matrix3::Identity();
  CHECK(foam::jacobian_determinant(F) == doctest::Approx(8.0));
}

TEST_CASE("green_lagrange_strain vanishes for F = I (no deformation)") {
  foam::Matrix3 F = foam::Matrix3::Identity();
  foam::Matrix3 E = foam::green_lagrange_strain(F);
  CHECK(E.norm() == doctest::Approx(0.0).epsilon(1e-12));
}

TEST_CASE("compute_deformation_gradient is a documented not-yet-implemented stub") {
  foam::MatrixX dummy(1, 1);
  dummy(0, 0) = 0.0;
  CHECK_THROWS_AS(foam::compute_deformation_gradient(dummy, dummy), std::logic_error);
}
