#include "doctest.h"
#include "foam/elements/quad4.hpp"
#include "foam/material_point/uniaxial_test.hpp"
#include "foam/materials/neo_hookean.hpp"
#include "foam/materials/ogden_hyperfoam.hpp"

namespace {

// Unit square element, nodes bottom-left/bottom-right/top-right/top-left
// (the convention Quad4Element and make_structured_grid_2d both use).
foam::MatrixX unit_square_positions() {
  foam::MatrixX ref_pos(4, 2);
  ref_pos << 0, 0, 1, 0, 1, 1, 0, 1;
  return ref_pos;
}

}  // namespace

TEST_CASE("Quad4Element: zero displacement gives zero internal force") {
  foam::NeoHookean material(1.0e5, 1.5e5);
  foam::Quad4Element element;
  foam::MatrixX disp = foam::MatrixX::Zero(4, 2);

  foam::VectorX f = element.internal_force(unit_square_positions(), disp, material);
  CHECK(f.norm() < 1e-8);
}

TEST_CASE("Quad4Element: homogeneous free-uniaxial-compression state matches the Step-1 "
          "single-material-point result exactly") {
  // This is the strongest available cross-check between the Step 1
  // (material-point) and Step 2 (FEM element) work: solve_uniaxial_point
  // finds the lateral stretch that makes lateral Cauchy stress vanish;
  // apply that SAME homogeneous F to a single unit-square element via
  // EXACT nodal displacements (u_a = (F-I)*X_a, no solving involved
  // here), and check that (a) the element's internal force has zero
  // lateral (y) component at every node -- since S22=0 exactly -- and
  // (b) the total x-force on the "loaded" edge equals the nominal axial
  // stress from Step 1 (both are force-per-original-area quantities).
  foam::OgdenHyperfoam material({foam::OgdenHyperfoam::Term{80000.0, 13.0, 0.0}});
  auto point = foam::solve_uniaxial_point(material, /*axial_stretch=*/0.7);
  REQUIRE(point.converged);

  foam::MatrixX ref_pos = unit_square_positions();
  foam::MatrixX disp(4, 2);
  for (int a = 0; a < 4; ++a) {
    disp(a, 0) = (point.axial_stretch - 1.0) * ref_pos(a, 0);
    disp(a, 1) = (point.lateral_stretch - 1.0) * ref_pos(a, 1);
  }

  foam::Quad4Element element;
  foam::VectorX f = element.internal_force(ref_pos, disp, material);

  // Lateral (y) components: indices 1,3,5,7 (node*2+1).
  for (int a = 0; a < 4; ++a) {
    CHECK(std::abs(f(2 * a + 1)) < 1e-6);
  }

  // Nodes 1 and 2 (bottom-right, top-right) sit on the x=1 edge.
  double right_edge_x_sum = f(2 * 1) + f(2 * 2);
  CHECK(right_edge_x_sum == doctest::Approx(point.nominal_stress_axial).epsilon(1e-6));
}

TEST_CASE("Quad4Element: tangent_stiffness matches a finite-difference derivative of "
          "internal_force (the element-level analogue of the material tangent check)") {
  foam::NeoHookean material(1.2e5, 1.8e5);
  foam::Quad4Element element;
  foam::MatrixX ref_pos = unit_square_positions();

  // A generic, moderately-deformed, non-symmetric displacement state
  // (not aligned with any special axis), so the check exercises the
  // general case rather than a special one.
  foam::MatrixX disp(4, 2);
  disp << 0.00, 0.00, -0.08, 0.02, -0.05, -0.12, 0.03, -0.07;

  foam::VectorX f0 = element.internal_force(ref_pos, disp, material);
  foam::MatrixX K = element.tangent_stiffness(ref_pos, disp, material);
  REQUIRE(K.rows() == 8);
  REQUIRE(K.cols() == 8);

  double h = 1e-6;
  for (int dof = 0; dof < 8; ++dof) {
    foam::MatrixX disp_plus = disp;
    foam::MatrixX disp_minus = disp;
    disp_plus(dof / 2, dof % 2) += h;
    disp_minus(dof / 2, dof % 2) -= h;

    foam::VectorX f_plus = element.internal_force(ref_pos, disp_plus, material);
    foam::VectorX f_minus = element.internal_force(ref_pos, disp_minus, material);
    foam::VectorX numerical_column = (f_plus - f_minus) / (2.0 * h);

    double scale = std::max(K.col(dof).cwiseAbs().maxCoeff(), 1.0);
    for (int row = 0; row < 8; ++row) {
      CHECK(std::abs(K(row, dof) - numerical_column(row)) < scale * 1e-3);
    }
  }
}

TEST_CASE("Quad4Element: rejects a degenerate/mis-ordered element (non-positive reference "
          "Jacobian)") {
  foam::NeoHookean material(1.0e5, 1.5e5);
  foam::Quad4Element element;
  foam::MatrixX bad_ref_pos(4, 2);
  // Same square, but listed CLOCKWISE instead of CCW -- this reverses
  // the sign of the Jacobian everywhere (guaranteed, unlike a
  // self-intersecting ordering which can have mixed-sign regions).
  bad_ref_pos << 0, 0, 0, 1, 1, 1, 1, 0;
  foam::MatrixX disp = foam::MatrixX::Zero(4, 2);

  CHECK_THROWS_AS(element.internal_force(bad_ref_pos, disp, material), std::runtime_error);
}
