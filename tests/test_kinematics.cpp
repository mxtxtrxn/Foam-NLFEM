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

TEST_CASE("compute_deformation_gradient: zero displacement gives F = I (2D and 3D)") {
  foam::MatrixX grad2(3, 2);
  grad2 << 0.5, 0.0, 0.0, 0.5, -0.5, -0.5;
  foam::MatrixX zero_disp2 = foam::MatrixX::Zero(3, 2);
  CHECK(foam::compute_deformation_gradient(grad2, zero_disp2).isApprox(foam::Matrix3::Identity()));

  foam::MatrixX grad3(4, 3);
  grad3 << 0.5, 0.0, 0.1, 0.0, 0.5, 0.1, -0.5, -0.5, 0.1, 0.0, 0.0, -0.3;
  foam::MatrixX zero_disp3 = foam::MatrixX::Zero(4, 3);
  CHECK(foam::compute_deformation_gradient(grad3, zero_disp3).isApprox(foam::Matrix3::Identity()));
}

TEST_CASE("compute_deformation_gradient: dim=2 matches a hand-computed example and embeds "
          "plane strain (F(2,2) = 1, rest of row/col 2 = 0)") {
  foam::MatrixX grad(3, 2);
  grad << 0.5, 0.0, 0.0, 0.5, -0.5, -0.5;
  foam::MatrixX disp(3, 2);
  disp << 0.1, 0.2, 0.3, 0.1, -0.05, 0.0;

  foam::Matrix3 F = foam::compute_deformation_gradient(grad, disp);

  // F(i,j) = delta_ij + sum_a disp(a,i)*grad(a,j), computed by hand:
  CHECK(F(0, 0) == doctest::Approx(1.075));
  CHECK(F(0, 1) == doctest::Approx(0.175));
  CHECK(F(1, 0) == doctest::Approx(0.1));
  CHECK(F(1, 1) == doctest::Approx(1.05));
  CHECK(F(2, 2) == doctest::Approx(1.0));
  CHECK(F(0, 2) == doctest::Approx(0.0));
  CHECK(F(1, 2) == doctest::Approx(0.0));
  CHECK(F(2, 0) == doctest::Approx(0.0));
  CHECK(F(2, 1) == doctest::Approx(0.0));
}

TEST_CASE("compute_deformation_gradient: rejects mismatched shape_gradients/nodal_displacements "
          "dimensions") {
  foam::MatrixX grad(3, 2);
  grad.setZero();
  foam::MatrixX disp(2, 2);  // wrong row count
  disp.setZero();
  CHECK_THROWS_AS(foam::compute_deformation_gradient(grad, disp), std::invalid_argument);
}
