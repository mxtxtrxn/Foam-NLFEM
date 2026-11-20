#pragma once

#include <functional>

#include "foam/core/types.hpp"

namespace foam {

// Generic Newton-Raphson iteration on an abstract nonlinear system
// R(x) = 0, with tangent J(x) = dR/dx. It knows nothing about meshes,
// elements, or materials -- Problem will wire it up to
// Assembler::assemble_residual/assemble_tangent later (roadmap step 2).
// Being decoupled like this means it can be implemented and unit-tested
// right now, against any toy nonlinear function, independent of whether
// the FEM/material layers above it are finished yet.
class NewtonRaphsonSolver {
 public:
  struct Options {
    Real tolerance = 1e-10;
    int max_iterations = 50;
    bool verbose = false;
  };

  struct Result {
    VectorX solution;
    bool converged = false;
    int iterations = 0;
    Real final_residual_norm = 0.0;
  };

  using ResidualFunction = std::function<VectorX(const VectorX&)>;
  using TangentFunction = std::function<SparseMatrix(const VectorX&)>;

  NewtonRaphsonSolver();
  explicit NewtonRaphsonSolver(Options options);

  // Solves R(x) = 0 starting from `initial_guess`, using `tangent` as the
  // Jacobian at each iterate. Does not throw on non-convergence; check
  // Result::converged instead (foam simulations routinely need to retry
  // with a smaller load step on non-convergence, not crash).
  Result solve(const ResidualFunction& residual, const TangentFunction& tangent,
               const VectorX& initial_guess) const;

 private:
  Options options_;
};

}  // namespace foam
