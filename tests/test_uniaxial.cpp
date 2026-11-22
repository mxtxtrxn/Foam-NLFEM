#include "doctest.h"
#include "foam/material_point/uniaxial_test.hpp"
#include "foam/materials/ogden_hyperfoam.hpp"

TEST_CASE("solve_uniaxial_point: beta=0 solves to lateral_stretch ~= 1 (zero Poisson's ratio)") {
  foam::OgdenHyperfoam material({foam::OgdenHyperfoam::Term{80000.0, 13.0, 0.0}});
  auto point = foam::solve_uniaxial_point(material, /*axial_stretch=*/0.7);

  CHECK(point.converged);
  CHECK(point.lateral_stretch == doctest::Approx(1.0).epsilon(1e-6));
  CHECK(point.cauchy_stress_axial < 0.0);  // compressive
}

TEST_CASE("solve_uniaxial_point: beta>0 bulges outward (lateral_stretch > 1) under compression") {
  foam::OgdenHyperfoam material({foam::OgdenHyperfoam::Term{80000.0, 13.0, 0.75}});  // nu = 0.3
  auto point = foam::solve_uniaxial_point(material, /*axial_stretch=*/0.7);

  CHECK(point.converged);
  CHECK(point.lateral_stretch > 1.0);
  CHECK(point.cauchy_stress_axial < 0.0);
}

TEST_CASE("uniaxial_compression_curve: produces a monotonically increasing compressive "
          "nominal-stress magnitude as strain increases") {
  foam::OgdenHyperfoam material({foam::OgdenHyperfoam::Term{80000.0, 13.0, 0.0}});
  auto curve = foam::uniaxial_compression_curve(material, /*axial_stretch_min=*/0.4,
                                                 /*axial_stretch_max=*/1.0, /*num_steps=*/20);

  REQUIRE(curve.size() == 21);
  CHECK(curve.front().engineering_strain == doctest::Approx(0.0).epsilon(1e-9));
  CHECK(curve.back().engineering_strain == doctest::Approx(0.6).epsilon(1e-9));

  for (std::size_t i = 1; i < curve.size(); ++i) {
    CHECK(curve[i].converged);
    // nominal_stress_axial is negative (compressive); its magnitude
    // should grow monotonically with strain for this material.
    CHECK(std::abs(curve[i].nominal_stress_axial) > std::abs(curve[i - 1].nominal_stress_axial));
  }
}

TEST_CASE("uniaxial_compression_curve: densification shows up as a LOCAL-stiffness upturn") {
  // NOTE: this deliberately does NOT use secant-from-origin stiffness
  // (stress/strain measured from e=0). For a toe/plateau/densification
  // curve, that secant naturally keeps DECREASING for a long stretch
  // after the toe -- the plateau's nearly-flat stress drags the running
  // average down even while the material is, locally, already starting
  // to densify again. Hand-deriving the closed-form beta=0 curve for
  // this exact material confirms it: secant stiffness at 80% strain
  // (~77 kPa) is still below its value at 5% strain (~126 kPa), even
  // though the material is unambiguously stiffer locally by then. The
  // correct diagnostic is LOCAL (point-to-point) stiffness, which is
  // what a real stress-strain plot's slope shows.
  foam::OgdenHyperfoam material({foam::OgdenHyperfoam::Term{80000.0, 13.0, 0.0}});
  auto curve = foam::uniaxial_compression_curve(material, 0.2, 1.0, 80);

  REQUIRE(curve.size() > 20);

  auto local_stiffness = [](const foam::UniaxialPoint& a, const foam::UniaxialPoint& b) {
    return std::abs(b.nominal_stress_axial - a.nominal_stress_axial) /
           std::abs(b.engineering_strain - a.engineering_strain);
  };

  double min_local_stiffness = 1e300;
  for (std::size_t i = 1; i < curve.size(); ++i) {
    min_local_stiffness = std::min(min_local_stiffness, local_stiffness(curve[i - 1], curve[i]));
  }
  double final_local_stiffness = local_stiffness(curve[curve.size() - 2], curve.back());

  CHECK(final_local_stiffness > 3.0 * min_local_stiffness);
}
