// Roadmap step 2: a real 2D mesh (Quad4 elements, Assembler,
// NewtonRaphsonSolver with Dirichlet BC elimination) under free
// uniaxial compression, cross-checked against Step 1's single
// -material-point result (foam::uniaxial_compression_curve). For a
// symmetric free-compression setup the deformation is homogeneous, so
// the two completely different code paths (one material point vs. a
// full mesh) should agree closely -- this is the strongest validation
// available for the whole assembled pipeline.

#include <cmath>
#include <iomanip>
#include <iostream>

#include "foam/elements/quad4.hpp"
#include "foam/io/csv_writer.hpp"
#include "foam/material_point/uniaxial_test.hpp"
#include "foam/materials/ogden_hyperfoam.hpp"
#include "foam/mesh/mesh_generators.hpp"
#include "foam/simulation/problem.hpp"

int main() {
  using foam::OgdenHyperfoam;

  // Same illustrative material as apps/single_element_test.cpp.
  OgdenHyperfoam material({OgdenHyperfoam::Term{80000.0, 13.0, 0.0}});

  foam::Real height = 1.0;
  foam::Real width = 1.0;
  int nx = 3, ny = 6;  // a somewhat fine mesh, to also show mesh-size independence
  foam::Real engineering_strain_target = 0.5;  // 50% compression
  foam::Real top_displacement = -engineering_strain_target * height;

  auto mesh = foam::make_structured_grid_2d(width, height, nx, ny);
  auto bcs = foam::uniaxial_compression_bcs(nx, ny, top_displacement);
  foam::Quad4Element element;
  foam::NeumannBCs no_neumann;

  std::cout << "Mesh: " << nx << "x" << ny << " Quad4 elements (" << mesh.num_nodes()
            << " nodes, " << mesh.num_elements() << " elements)\n";
  std::cout << "Free uniaxial compression to " << (engineering_strain_target * 100.0)
            << "% engineering strain\n\n";

  foam::Problem problem(mesh, material, element, bcs, no_neumann);
  auto fem_curve = problem.run_compression_test(/*num_steps=*/40);

  if (fem_curve.size() < 41) {
    std::cout << "NOTE: the load-stepping sweep stopped early at "
              << (fem_curve.back().first * 100.0)
              << "% strain (Newton did not converge past this point).\n\n";
  }

  // Step 1 reference: same material, same axial-stretch history, a
  // single material point instead of a mesh.
  auto point_curve = foam::uniaxial_compression_curve(
      material, /*axial_stretch_min=*/1.0 - engineering_strain_target, /*axial_stretch_max=*/1.0,
      /*num_steps=*/40);

  std::cout << std::fixed << std::setprecision(2);
  std::cout << std::setw(10) << "strain(%)" << std::setw(18) << "FEM stress(kPa)" << std::setw(22)
            << "material-point(kPa)" << std::setw(14) << "diff(%)\n";

  std::vector<std::pair<foam::Real, foam::Real>> fem_csv;
  for (std::size_t i = 0; i < fem_curve.size() && i < point_curve.size(); ++i) {
    foam::Real fem_stress = fem_curve[i].second;
    foam::Real point_stress = point_curve[i].nominal_stress_axial;
    foam::Real diff_pct =
        std::abs(point_stress) > 1.0 ? 100.0 * std::abs(fem_stress - point_stress) /
                                            std::abs(point_stress)
                                     : 0.0;
    if (i % 5 == 0 || i == fem_curve.size() - 1) {
      std::cout << std::setw(10) << (fem_curve[i].first * 100.0) << std::setw(18)
                << (fem_stress / 1000.0) << std::setw(22) << (point_stress / 1000.0)
                << std::setw(13) << diff_pct << "%\n";
    }
    fem_csv.emplace_back(fem_curve[i].first, fem_stress);
  }

  foam::write_csv("mesh_compression_curve.csv", fem_csv, "engineering_strain", "nominal_stress_pa");
  std::cout << "\n-> curve written to mesh_compression_curve.csv\n";

  return 0;
}
