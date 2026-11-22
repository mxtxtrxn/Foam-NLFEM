#include "doctest.h"
#include "foam/elements/quad4.hpp"
#include "foam/material_point/uniaxial_test.hpp"
#include "foam/materials/ogden_hyperfoam.hpp"
#include "foam/mesh/mesh_generators.hpp"
#include "foam/simulation/problem.hpp"

TEST_CASE("uniaxial_compression_bcs: bottom row y=0, top row y=target, one x pinned") {
  auto bcs = foam::uniaxial_compression_bcs(/*nx=*/2, /*ny=*/3, /*top_y_displacement=*/-0.4);

  int bottom_count = 0, top_count = 0, x_pins = 0;
  for (const auto& bc : bcs) {
    if (bc.dof == 1 && bc.value == 0.0) ++bottom_count;
    if (bc.dof == 1 && bc.value == -0.4) ++top_count;
    if (bc.dof == 0) ++x_pins;
  }
  CHECK(bottom_count == 3);  // nx+1 bottom nodes
  CHECK(top_count == 3);     // nx+1 top nodes
  CHECK(x_pins == 1);
}

TEST_CASE("Problem::run_compression_test: a full 2D mesh under free uniaxial compression "
          "matches Step 1's single-material-point curve (the end-to-end cross-check)") {
  // Same material, same axial strain history, two completely different
  // code paths: Step 1 solves one material point directly; here, a real
  // mesh goes through Assembler + Quad4Element + NewtonRaphsonSolver
  // with Dirichlet BC elimination. For a homogeneous deformation state
  // (which this free-compression setup produces, by symmetry) they
  // must agree, up to Newton/discretization tolerance.
  foam::OgdenHyperfoam material({foam::OgdenHyperfoam::Term{80000.0, 13.0, 0.0}});

  foam::Real height = 1.0;
  int nx = 2, ny = 3;
  foam::Real engineering_strain_target = 0.3;  // 30% compression
  foam::Real top_displacement = -engineering_strain_target * height;

  auto mesh = foam::make_structured_grid_2d(1.0, height, nx, ny);
  auto bcs = foam::uniaxial_compression_bcs(nx, ny, top_displacement);
  foam::Quad4Element element;
  foam::NeumannBCs no_neumann;

  foam::Problem problem(mesh, material, element, bcs, no_neumann);
  auto fem_curve = problem.run_compression_test(/*num_steps=*/10);

  REQUIRE(fem_curve.size() == 11);
  CHECK(fem_curve.back().first == doctest::Approx(engineering_strain_target).epsilon(1e-6));

  // Cross-check the final point against the single-material-point result.
  auto point = foam::solve_uniaxial_point(material, 1.0 - engineering_strain_target);
  REQUIRE(point.converged);

  CHECK(fem_curve.back().second ==
        doctest::Approx(point.nominal_stress_axial).epsilon(1e-4));
}
