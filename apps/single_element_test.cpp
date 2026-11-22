// Roadmap step 1: validate OgdenHyperfoam on a single material point,
// under the boundary condition that actually matches a real compression
// test -- free lateral surfaces (unconfined/"free" uniaxial compression),
// not lateral stretch pinned to 1 (that would be uniaxial STRAIN, a
// confined and much easier-to-get-right-by-accident test case).
//
// IMPORTANT: the material parameters below are ILLUSTRATIVE. They were
// chosen to produce the qualitatively correct three-regime foam
// response (elastic toe / plateau / densification) discussed in
// docs/roadmap.md, not fitted to one specific cited experimental
// dataset. Fitting against real test data is roadmap step 4.

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "foam/io/csv_writer.hpp"
#include "foam/material_point/uniaxial_test.hpp"
#include "foam/materials/ogden_hyperfoam.hpp"

namespace {

void print_header(const std::string& title) {
  std::cout << "\n=== " << title << " ===\n";
}

void report_curve(const std::string& label, const std::vector<foam::UniaxialPoint>& curve,
                   const std::string& csv_path) {
  print_header(label);

  if (curve.size() < 2) {
    std::cout << "  (too few converged points to report)\n";
    return;
  }
  if (!curve.back().converged) {
    std::cout << "  NOTE: sweep stopped early at " << (curve.back().engineering_strain * 100.0)
               << "% engineering strain (Newton did not converge past this point --\n"
               << "  expected with a plain full-Newton step deep in densification; see\n"
               << "  docs/roadmap.md step 3, arc-length / step limiting).\n";
  }

  // Initial tangent modulus, from the first two points (small-strain).
  const auto& p0 = curve[0];
  const auto& p1 = curve[1];
  double initial_modulus_kpa =
      (p1.nominal_stress_axial - p0.nominal_stress_axial) /
      (p1.engineering_strain - p0.engineering_strain) / 1000.0;

  std::cout << std::fixed << std::setprecision(3);
  std::cout << "  initial tangent modulus (elastic toe) : " << initial_modulus_kpa << " kPa\n";

  auto report_at = [&](double target_strain) {
    const foam::UniaxialPoint* closest = nullptr;
    double best_diff = 1e9;
    for (const auto& p : curve) {
      double diff = std::abs(p.engineering_strain - target_strain);
      if (diff < best_diff) {
        best_diff = diff;
        closest = &p;
      }
    }
    if (closest) {
      std::cout << "  @ " << std::setw(5) << (closest->engineering_strain * 100.0)
                << "% strain: nominal stress = " << std::setw(9)
                << (closest->nominal_stress_axial / 1000.0)
                << " kPa, lateral stretch = " << closest->lateral_stretch << "\n";
    }
  };

  report_at(0.10);
  report_at(0.40);
  report_at(curve.back().engineering_strain);

  foam::write_csv(csv_path, [&] {
    std::vector<std::pair<foam::Real, foam::Real>> xy;
    xy.reserve(curve.size());
    for (const auto& p : curve) {
      xy.emplace_back(p.engineering_strain, p.nominal_stress_axial);
    }
    return xy;
  }(),
                  "engineering_strain", "nominal_stress_pa");
  std::cout << "  -> curve written to " << csv_path << "\n";
}

}  // namespace

int main() {
  using foam::OgdenHyperfoam;

  // Single-term (N=1) Ogden hyperfoam, two variants differing only in
  // beta (equivalently, in the effective Poisson's ratio):
  //   beta = 0     -> zero Poisson's ratio (common default assumption
  //                   for many open-cell foams; see e.g. Abaqus foam
  //                   -fitting examples, which often fix beta_i = 0)
  //   beta = 0.75  -> nu = beta/(1+2*beta) ... actually nu solves
  //                   beta = nu/(1-2*nu), so nu = beta/(1+2*beta) = 0.3,
  //                   representative of a foam with noticeable lateral
  //                   bulging under compression (closed-cell-like).
  OgdenHyperfoam foam_zero_poisson({OgdenHyperfoam::Term{80000.0, 13.0, 0.0}});
  OgdenHyperfoam foam_positive_poisson({OgdenHyperfoam::Term{80000.0, 13.0, 0.75}});

  print_header("Sanity check: undeformed state gives zero stress");
  foam::Matrix3 I = foam::Matrix3::Identity();
  std::cout << "  ||sigma(F=I)|| (beta=0)    = " << foam_zero_poisson.cauchy_stress(I).norm()
            << "\n";
  std::cout << "  ||sigma(F=I)|| (beta=0.75) = " << foam_positive_poisson.cauchy_stress(I).norm()
            << "\n";

  // Sweep from the undeformed state (axial_stretch = 1) down to 70%
  // engineering compressive strain (axial_stretch = 0.3).
  auto curve_zero_poisson =
      foam::uniaxial_compression_curve(foam_zero_poisson, /*axial_stretch_min=*/0.3,
                                        /*axial_stretch_max=*/1.0, /*num_steps=*/60);
  auto curve_positive_poisson =
      foam::uniaxial_compression_curve(foam_positive_poisson, 0.3, 1.0, 60);

  report_curve("beta = 0 (zero Poisson's ratio)", curve_zero_poisson,
               "foam_beta0_curve.csv");
  report_curve("beta = 0.75 (nu = 0.3)", curve_positive_poisson,
               "foam_beta_pos_curve.csv");

  std::cout << "\nNote on lateral behavior: under axial COMPRESSION, beta=0 keeps the\n"
               "lateral stretch essentially at 1 (no bulging), while beta=0.75 shows\n"
               "lateral stretch > 1 (the specimen bulges outward), as expected for a\n"
               "material with a positive effective Poisson's ratio.\n";

  return 0;
}
