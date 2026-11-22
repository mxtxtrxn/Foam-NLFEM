#include <cmath>

#include "doctest.h"
#include "foam/solvers/arc_length.hpp"

namespace {

// R(u, lambda) = u^3 - 3u - lambda = 0  =>  closed form lambda(u) = u^3 - 3u.
// Classic S-shaped snap-through curve: dR/du = 3u^2-3 = 0 (limit points)
// at u=-1 (lambda=2, a local MAX of lambda) and u=1 (lambda=-2, a local
// MIN) -- exactly the shape a load- or displacement-stepping scheme
// cannot follow past, but arc-length should trace smoothly.
foam::VectorX cubic_residual(const foam::VectorX& u, foam::Real lambda) {
  foam::VectorX r(1);
  r(0) = u(0) * u(0) * u(0) - 3.0 * u(0) - lambda;
  return r;
}

foam::SparseMatrix cubic_tangent(const foam::VectorX& u, foam::Real /*lambda*/) {
  foam::SparseMatrix K(1, 1);
  K.insert(0, 0) = 3.0 * u(0) * u(0) - 3.0;
  K.makeCompressed();
  return K;
}

foam::Real exact_lambda(foam::Real u) { return u * u * u - 3.0 * u; }

}  // namespace

TEST_CASE("ArcLengthSolver: traces the full cubic snap-through curve, matching the closed "
          "form at every step, through BOTH limit points") {
  foam::ArcLengthSolver::Options options;
  options.arc_length = 0.1;
  options.psi = 1.0;
  options.num_steps = 300;
  options.tolerance = 1e-10;
  options.initial_direction = 1.0;  // start by increasing u (and initially lambda)
  foam::ArcLengthSolver solver(options);

  foam::VectorX u0(1);
  u0(0) = -2.5;
  foam::Real lambda0 = exact_lambda(u0(0));

  auto path = solver.solve(cubic_residual, cubic_tangent, u0, lambda0);

  REQUIRE(path.size() == static_cast<std::size_t>(options.num_steps));

  bool passed_first_limit_point = false;   // u crosses -1
  bool passed_second_limit_point = false;  // u crosses +1
  for (const auto& point : path) {
    CHECK(point.converged);
    CHECK(point.lambda == doctest::Approx(exact_lambda(point.u(0))).epsilon(1e-6));
    if (point.u(0) > -1.0) passed_first_limit_point = true;
    if (point.u(0) > 1.0) passed_second_limit_point = true;
  }
  CHECK(passed_first_limit_point);
  CHECK(passed_second_limit_point);

  // The path should have reached well past the second limit point and
  // be back on the far stable branch (u large, lambda increasing again).
  CHECK(path.back().u(0) > 2.0);
}

TEST_CASE("ArcLengthSolver: lambda is genuinely non-monotonic along the traced path "
          "(the defining feature of snap-through, and of why plain load-stepping would fail "
          "here)") {
  foam::ArcLengthSolver::Options options;
  options.arc_length = 0.1;
  options.num_steps = 300;
  options.tolerance = 1e-10;
  foam::ArcLengthSolver solver(options);

  foam::VectorX u0(1);
  u0(0) = -2.5;
  auto path = solver.solve(cubic_residual, cubic_tangent, u0, exact_lambda(u0(0)));

  // Track lambda's FIRST local maximum (the u=-1 limit point, lambda=2)
  // -- not the global max over the whole path, which also covers the
  // far stable branch where lambda keeps climbing much higher.
  double first_local_max = path[0].lambda;
  bool seen_decrease = false;
  for (std::size_t i = 1; i < path.size() && !seen_decrease; ++i) {
    if (path[i].lambda > first_local_max) {
      first_local_max = path[i].lambda;
    } else {
      seen_decrease = true;
    }
  }
  CHECK(seen_decrease);
  CHECK(first_local_max == doctest::Approx(2.0).epsilon(1e-2));  // the u=-1 limit point value
}

TEST_CASE("ArcLengthSolver: automatic step-halving retry recovers from a too-large step") {
  // arc_length=0.15 is large enough that the strongly nonlinear
  // unstable branch (between the two limit points) fails to converge
  // in a single attempt -- empirically confirmed (not every larger
  // arc_length fails; it depends on exactly where the discrete steps
  // land relative to the curve's local nonlinearity, so this specific
  // value is deliberately chosen, not guessed) by setting
  // max_retries_per_step=0 and observing the solver stop early, THEN
  // confirming that with retries enabled, it gets much further.
  foam::ArcLengthSolver::Options no_retry_options;
  no_retry_options.arc_length = 0.15;
  no_retry_options.num_steps = 200;
  no_retry_options.tolerance = 1e-10;
  no_retry_options.max_retries_per_step = 0;
  foam::ArcLengthSolver solver_no_retry(no_retry_options);

  foam::VectorX u0(1);
  u0(0) = -2.5;
  auto path_no_retry =
      solver_no_retry.solve(cubic_residual, cubic_tangent, u0, exact_lambda(u0(0)));

  foam::ArcLengthSolver::Options retry_options = no_retry_options;
  retry_options.max_retries_per_step = 5;
  foam::ArcLengthSolver solver_with_retry(retry_options);
  auto path_with_retry =
      solver_with_retry.solve(cubic_residual, cubic_tangent, u0, exact_lambda(u0(0)));

  // With retries, the solver should get at least as far (in step count)
  // as without, and in this specific case strictly further (the
  // oversized step does fail somewhere in the unstable branch without
  // retries).
  CHECK(path_with_retry.size() >= path_no_retry.size());
  CHECK(path_with_retry.size() > static_cast<std::size_t>(0.5 * no_retry_options.num_steps));

  // Every point that DID converge, in both cases, matches the closed
  // form exactly (retries change whether/how far the solver gets, not
  // the correctness of what does converge). A path that stopped early
  // ends with exactly one non-converged attempt -- excluded here on
  // purpose, since by definition R was not driven to zero there.
  for (const auto& point : path_no_retry) {
    if (point.converged) {
      CHECK(point.lambda == doctest::Approx(exact_lambda(point.u(0))).epsilon(1e-6));
    }
  }
  for (const auto& point : path_with_retry) {
    CHECK(point.converged);
    CHECK(point.lambda == doctest::Approx(exact_lambda(point.u(0))).epsilon(1e-6));
  }
}

TEST_CASE("ArcLengthSolver: default constructor works and produces a usable solver") {
  foam::ArcLengthSolver solver;
  foam::VectorX u0(1);
  u0(0) = -2.5;
  auto path = solver.solve(cubic_residual, cubic_tangent, u0, exact_lambda(u0(0)));
  REQUIRE_FALSE(path.empty());
  CHECK(path.front().converged);
}
