#include <cmath>

#include "doctest.h"
#include "foam/solvers/linear_solver.hpp"
#include "foam/solvers/newton_raphson.hpp"

TEST_CASE("LinearSolver: solves a small diagonal SPD system exactly") {
  foam::SparseMatrix K(2, 2);
  K.insert(0, 0) = 2.0;
  K.insert(1, 1) = 3.0;
  K.makeCompressed();

  foam::VectorX rhs(2);
  rhs << 4.0, 9.0;

  foam::VectorX x = foam::LinearSolver::solve(K, rhs);
  CHECK(x(0) == doctest::Approx(2.0));
  CHECK(x(1) == doctest::Approx(3.0));
}

TEST_CASE("NewtonRaphsonSolver: converges on the scalar equation x^2 = 2") {
  foam::NewtonRaphsonSolver::Options options;
  options.tolerance = 1e-12;
  foam::NewtonRaphsonSolver solver(options);

  auto residual = [](const foam::VectorX& x) -> foam::VectorX {
    foam::VectorX r(1);
    r(0) = x(0) * x(0) - 2.0;
    return r;
  };
  auto tangent = [](const foam::VectorX& x) -> foam::SparseMatrix {
    foam::SparseMatrix J(1, 1);
    J.insert(0, 0) = 2.0 * x(0);
    J.makeCompressed();
    return J;
  };

  foam::VectorX initial_guess(1);
  initial_guess(0) = 1.0;

  auto result = solver.solve(residual, tangent, initial_guess);
  CHECK(result.converged);
  CHECK(result.solution(0) == doctest::Approx(std::sqrt(2.0)));
}

TEST_CASE("NewtonRaphsonSolver: reports non-convergence instead of throwing") {
  foam::NewtonRaphsonSolver::Options options;
  options.max_iterations = 1;  // deliberately too few to converge
  options.tolerance = 1e-14;
  foam::NewtonRaphsonSolver solver(options);

  auto residual = [](const foam::VectorX& x) -> foam::VectorX {
    foam::VectorX r(1);
    r(0) = x(0) * x(0) - 2.0;
    return r;
  };
  auto tangent = [](const foam::VectorX& x) -> foam::SparseMatrix {
    foam::SparseMatrix J(1, 1);
    J.insert(0, 0) = 2.0 * x(0);
    J.makeCompressed();
    return J;
  };

  foam::VectorX initial_guess(1);
  initial_guess(0) = 1000.0;  // far from the root, 1 iteration won't be enough

  auto result = solver.solve(residual, tangent, initial_guess);
  CHECK_FALSE(result.converged);
}
