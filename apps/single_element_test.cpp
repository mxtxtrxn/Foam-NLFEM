// NOTE: This is a placeholder, not the real roadmap-step-1 test.
//
// It solves the toy scalar equation x^2 = 2 with the generic
// NewtonRaphsonSolver, to prove the solver + sparse linear algebra
// plumbing works end to end *before* any foam material physics exists.
//
// Once OgdenHyperfoam::cauchy_stress/tangent_moduli are implemented,
// this file is meant to be replaced by the real single-element
// homogeneous-compression test described in the roadmap: build
// make_single_hex_element(), prescribe a compressive displacement on the
// top face, and compare the resulting nominal stress-strain curve
// against a published/analytical Ogden hyperfoam curve.

#include <cmath>
#include <iostream>

#include "foam/solvers/newton_raphson.hpp"

int main() {
  using foam::NewtonRaphsonSolver;
  using foam::SparseMatrix;
  using foam::VectorX;

  NewtonRaphsonSolver::Options options;
  options.verbose = true;
  NewtonRaphsonSolver solver(options);

  auto residual = [](const VectorX& x) -> VectorX {
    VectorX r(1);
    r(0) = x(0) * x(0) - 2.0;
    return r;
  };

  auto tangent = [](const VectorX& x) -> SparseMatrix {
    SparseMatrix J(1, 1);
    J.insert(0, 0) = 2.0 * x(0);
    J.makeCompressed();
    return J;
  };

  VectorX initial_guess(1);
  initial_guess(0) = 1.0;

  std::cout << "Sanity check: solving x^2 = 2 with the generic Newton-Raphson solver\n";
  std::cout << "(this only exercises the solver infrastructure, not foam physics yet)\n\n";

  auto result = solver.solve(residual, tangent, initial_guess);

  std::cout << "\nx            = " << result.solution(0) << "\n";
  std::cout << "expected     = " << std::sqrt(2.0) << "\n";
  std::cout << "converged    = " << (result.converged ? "yes" : "no") << "\n";
  std::cout << "iterations   = " << result.iterations << "\n";

  return result.converged ? 0 : 1;
}
