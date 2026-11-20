#include "foam/solvers/newton_raphson.hpp"

#include <cstdio>

#include "foam/solvers/linear_solver.hpp"

namespace foam {

NewtonRaphsonSolver::NewtonRaphsonSolver() : NewtonRaphsonSolver(Options{}) {}

NewtonRaphsonSolver::NewtonRaphsonSolver(Options options) : options_(options) {}

NewtonRaphsonSolver::Result NewtonRaphsonSolver::solve(const ResidualFunction& residual,
                                                        const TangentFunction& tangent,
                                                        const VectorX& initial_guess) const {
  Result result;
  result.solution = initial_guess;

  for (int iter = 0; iter < options_.max_iterations; ++iter) {
    VectorX R = residual(result.solution);
    Real residual_norm = R.norm();
    result.iterations = iter;
    result.final_residual_norm = residual_norm;

    if (options_.verbose) {
      std::printf("  Newton iter %2d: ||R|| = %.6e\n", iter, residual_norm);
    }

    if (residual_norm < options_.tolerance) {
      result.converged = true;
      return result;
    }

    SparseMatrix J = tangent(result.solution);
    VectorX delta = LinearSolver::solve(J, -R);
    result.solution += delta;
  }

  result.converged = false;
  return result;
}

}  // namespace foam
