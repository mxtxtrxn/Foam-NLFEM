// Roadmap step 3: demonstrates ArcLengthSolver on the classic
// textbook snap-through example, R(u,lambda) = u^3 - 3u - lambda = 0,
// whose closed form lambda(u) = u^3 - 3u is an S-shaped curve with
// limit points at u=-1 (lambda=2) and u=+1 (lambda=-2).
//
// IMPORTANT, and worth being upfront about: this is a toy 1-DOF
// problem, not the foam model. Checking apps/single_element_test's
// material (Ogden hyperfoam, beta=0, mu=80kPa, alpha=13) shows its
// uniaxial compression curve has no actual snap-through -- local
// stiffness drops through the plateau but never goes negative (see
// docs/roadmap.md). So plain displacement-controlled Newton-Raphson
// (already working end to end since roadmap step 2) is enough for
// THAT material; this demo validates the solver itself, against a
// problem where the difference between plain Newton and arc-length is
// unambiguous and checkable against a closed-form answer. Applying
// arc-length to an actually-unstable foam scenario (a softening
// material law, or geometric buckling in a finer mesh) is a natural
// next step, not done here.

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "foam/io/csv_writer.hpp"
#include "foam/solvers/arc_length.hpp"
#include "foam/solvers/newton_raphson.hpp"

namespace {

foam::VectorX residual(const foam::VectorX& u, foam::Real lambda) {
  foam::VectorX r(1);
  r(0) = u(0) * u(0) * u(0) - 3.0 * u(0) - lambda;
  return r;
}

foam::SparseMatrix tangent(const foam::VectorX& u, foam::Real /*lambda*/) {
  foam::SparseMatrix K(1, 1);
  K.insert(0, 0) = 3.0 * u(0) * u(0) - 3.0;
  K.makeCompressed();
  return K;
}

foam::Real exact_lambda(foam::Real u) { return u * u * u - 3.0 * u; }

}  // namespace

int main() {
  foam::VectorX u0(1);
  u0(0) = -2.5;
  foam::Real lambda0 = exact_lambda(u0(0));

  // --- Plain load control: increment lambda, solve for u each time ---
  std::cout << "=== Plain load-controlled Newton-Raphson ===\n";
  std::cout << "Stepping lambda from " << lambda0 << " upward in fixed increments...\n";
  foam::Real u_current = u0(0);
  foam::Real lambda_step = 0.2;
  int load_control_steps_survived = 0;
  for (int i = 1; i <= 60; ++i) {
    foam::Real lambda_target = lambda0 + i * lambda_step;
    foam::NewtonRaphsonSolver::Options options;
    options.tolerance = 1e-10;
    options.max_iterations = 50;
    foam::NewtonRaphsonSolver solver(options);

    foam::VectorX guess(1);
    guess(0) = u_current;
    auto result = solver.solve(
        [lambda_target](const foam::VectorX& u) { return residual(u, lambda_target); },
        [](const foam::VectorX& u) { return tangent(u, 0.0); }, guess);

    if (!result.converged) {
      std::cout << "  FAILED to converge at lambda=" << lambda_target
                << " (started from u=" << u_current << ") -- this is the limit point at "
                << "lambda=2, u=-1 that load control cannot pass.\n";
      break;
    }
    u_current = result.solution(0);
    load_control_steps_survived = i;
  }
  std::cout << "Plain load control survived " << load_control_steps_survived
            << " steps before failing.\n\n";

  // --- Arc-length: traces the FULL S-curve, including both limit points ---
  std::cout << "=== ArcLengthSolver ===\n";
  foam::ArcLengthSolver::Options arc_options;
  arc_options.arc_length = 0.1;
  arc_options.num_steps = 220;
  arc_options.tolerance = 1e-10;
  foam::ArcLengthSolver arc_solver(arc_options);

  auto path = arc_solver.solve(residual, tangent, u0, lambda0);

  std::cout << "Traced " << path.size() << " / " << arc_options.num_steps << " steps.\n";
  std::cout << std::fixed << std::setprecision(4);
  std::cout << std::setw(8) << "step" << std::setw(10) << "u" << std::setw(10) << "lambda"
            << std::setw(14) << "exact_lambda\n";

  std::vector<std::pair<foam::Real, foam::Real>> csv_data;
  for (std::size_t i = 0; i < path.size(); ++i) {
    foam::Real u = path[i].u(0);
    foam::Real exact = exact_lambda(u);
    if (i % 15 == 0 || i == path.size() - 1) {
      std::cout << std::setw(8) << i << std::setw(10) << u << std::setw(10) << path[i].lambda
                << std::setw(14) << exact << "\n";
    }
    csv_data.emplace_back(u, path[i].lambda);
  }

  bool passed_both_limit_points =
      path.back().u(0) > 1.0 && std::any_of(path.begin(), path.end(), [](const auto& p) {
        return p.u(0) < -1.0;
      });
  std::cout << "\nPassed both limit points (u=-1 and u=+1): "
            << (passed_both_limit_points ? "yes" : "no") << "\n";

  foam::write_csv("snap_through_curve.csv", csv_data, "u", "lambda");
  std::cout << "-> curve written to snap_through_curve.csv\n";

  return 0;
}
