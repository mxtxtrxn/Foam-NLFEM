#include "foam/material_point/uniaxial_test.hpp"

#include <algorithm>
#include <cmath>

#include "foam/solvers/newton_raphson.hpp"

namespace foam {

namespace {

Real lateral_stress(const MaterialModel& material, Real axial, Real lateral) {
  Matrix3 F = Matrix3::Zero();
  F(0, 0) = axial;
  F(1, 1) = lateral;
  F(2, 2) = lateral;
  return material.cauchy_stress(F)(1, 1);
}

}  // namespace

UniaxialPoint solve_uniaxial_point(const MaterialModel& material, Real axial_stretch,
                                    Real initial_lateral_guess) {
  // NewtonRaphsonSolver's stopping criterion is an ABSOLUTE residual
  // norm, but the residual here is itself a stress (Pa), and foam
  // moduli range from ~10 kPa to several MPa depending on the material
  // -- so scale the tolerance to the material's own stress level
  // instead of hard-coding a number that would be meaningless for some
  // materials and unreachable (due to floating-point noise) for others.
  Real reference_stress = std::abs(lateral_stress(material, axial_stretch, initial_lateral_guess));
  Real tolerance = std::max(reference_stress, Real(1.0)) * 1e-10;

  NewtonRaphsonSolver::Options options;
  options.tolerance = tolerance;
  options.max_iterations = 60;
  NewtonRaphsonSolver solver(options);

  auto residual = [&material, axial_stretch](const VectorX& x) -> VectorX {
    VectorX r(1);
    r(0) = lateral_stress(material, axial_stretch, x(0));
    return r;
  };

  auto tangent = [&material, axial_stretch](const VectorX& x) -> SparseMatrix {
    Real lateral = x(0);
    Real h = 1e-6 * std::max(std::abs(lateral), Real(1.0));
    Real r_plus = lateral_stress(material, axial_stretch, lateral + h);
    Real r_minus = lateral_stress(material, axial_stretch, lateral - h);
    SparseMatrix J(1, 1);
    J.insert(0, 0) = (r_plus - r_minus) / (2.0 * h);
    J.makeCompressed();
    return J;
  };

  VectorX x0(1);
  x0(0) = initial_lateral_guess;

  auto result = solver.solve(residual, tangent, x0);

  UniaxialPoint point;
  point.axial_stretch = axial_stretch;
  point.lateral_stretch = result.solution(0);
  point.engineering_strain = 1.0 - axial_stretch;
  point.converged = result.converged;

  Matrix3 F = Matrix3::Zero();
  F(0, 0) = axial_stretch;
  F(1, 1) = point.lateral_stretch;
  F(2, 2) = point.lateral_stretch;
  Matrix3 sigma = material.cauchy_stress(F);
  Real J = F.determinant();

  point.cauchy_stress_axial = sigma(0, 0);
  // Nominal (1st Piola-Kirchhoff) axial stress: P = J*sigma*F^{-T}, so
  // for this diagonal F, P_11 = J*sigma_11/axial_stretch -- this is
  // what a real test machine reports (force / ORIGINAL cross-section).
  point.nominal_stress_axial = J * sigma(0, 0) / axial_stretch;

  return point;
}

std::vector<UniaxialPoint> uniaxial_compression_curve(const MaterialModel& material,
                                                        Real axial_stretch_min,
                                                        Real axial_stretch_max, int num_steps) {
  std::vector<UniaxialPoint> curve;
  curve.reserve(static_cast<std::size_t>(num_steps) + 1);

  Real lateral_guess = 1.0;
  for (int i = 0; i <= num_steps; ++i) {
    Real t = static_cast<Real>(i) / static_cast<Real>(num_steps);
    Real axial = axial_stretch_max + t * (axial_stretch_min - axial_stretch_max);

    try {
      UniaxialPoint point = solve_uniaxial_point(material, axial, lateral_guess);
      lateral_guess = point.lateral_stretch;  // warm-start next step (continuation)
      curve.push_back(point);
      if (!point.converged) {
        break;  // report what converged rather than silently extrapolating
      }
    } catch (const std::exception&) {
      // Element inverted (J <= 0) or similar: the plain full-Newton
      // step in solve_uniaxial_point has no line search, so it can
      // overshoot into an invalid lateral stretch at large compressive
      // strains. Stop the sweep here rather than crash; this is
      // exactly the kind of failure arc-length/step-limiting
      // (roadmap step 3) is meant to fix.
      break;
    }
  }
  return curve;
}

}  // namespace foam
