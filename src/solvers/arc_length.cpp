#include "foam/solvers/arc_length.hpp"

#include <cmath>
#include <stdexcept>

#include "foam/solvers/linear_solver.hpp"

namespace foam {

namespace {

// Generic dR/dlambda via central finite difference -- same philosophy
// as numerical_tangent_moduli(), just one scalar direction instead of
// six.
VectorX numerical_dR_dlambda(const ArcLengthSolver::ResidualFunction& residual, const VectorX& u,
                              Real lambda, Real h) {
  VectorX r_plus = residual(u, lambda + h);
  VectorX r_minus = residual(u, lambda - h);
  return (r_plus - r_minus) / (2.0 * h);
}

// One arc-length predictor-corrector attempt at a GIVEN arc length
// (the caller retries with a smaller arc_length on failure). Returns
// the step's (delta_u, delta_lambda, converged, iterations); u_current
// /lambda_current are the already-converged state this step starts
// from, and prev_delta_u/prev_delta_lambda feed the sign-continuity
// criterion.
struct AttemptResult {
  VectorX delta_u;
  Real delta_lambda = 0.0;
  bool converged = false;
  int iterations = 0;
};

AttemptResult attempt_step(const ArcLengthSolver::ResidualFunction& residual,
                            const ArcLengthSolver::TangentFunction& tangent,
                            const ArcLengthSolver::Options& options, const VectorX& u_current,
                            Real lambda_current, const VectorX& prev_delta_u,
                            Real prev_delta_lambda, Real arc_length) {
  // --- Predictor: move along the tangent to the equilibrium path ---
  SparseMatrix K0 = tangent(u_current, lambda_current);
  VectorX q0 =
      numerical_dR_dlambda(residual, u_current, lambda_current, options.lambda_perturbation);
  // Along the path, d(R)/dlambda_total = K*(du/dlambda) + q = 0, so the
  // tangent direction is v0 = du/dlambda = -K^-1*q.
  VectorX v0 = LinearSolver::solve(K0, -q0);

  Real denom0 = std::sqrt(v0.squaredNorm() + options.psi * options.psi);
  if (denom0 < 1e-13) {
    throw std::runtime_error(
        "ArcLengthSolver: degenerate predictor direction (||v0|| and psi both ~0)");
  }

  // Sign: continue in the direction consistent with the previous
  // step's increment (automatically flips when passing a limit point,
  // since v0 itself changes character there).
  Real continuity = v0.dot(prev_delta_u) + options.psi * options.psi * prev_delta_lambda;
  Real sign = (continuity >= 0.0) ? 1.0 : -1.0;

  // Chosen so ||delta_u||^2 + psi^2*delta_lambda^2 = arc_length^2 with
  // delta_u = delta_lambda_predictor * v0.
  Real delta_lambda_predictor = sign * arc_length / denom0;
  VectorX delta_u = delta_lambda_predictor * v0;
  Real delta_lambda = delta_lambda_predictor;

  // --- Corrector: Newton iterations constrained to the hyperplane
  // through the predictor point, normal to (v0, psi^2*delta_lambda0)
  // -- the classical Riks/Wempner "normal plane" method. ---
  AttemptResult result;
  int iter = 0;
  for (; iter < options.max_iterations_per_step; ++iter) {
    VectorX u_trial = u_current + delta_u;
    Real lambda_trial = lambda_current + delta_lambda;
    VectorX R = residual(u_trial, lambda_trial);
    if (R.norm() < options.tolerance) {
      result.converged = true;
      break;
    }

    SparseMatrix K = tangent(u_trial, lambda_trial);
    VectorX q = numerical_dR_dlambda(residual, u_trial, lambda_trial, options.lambda_perturbation);

    VectorX du_bar1 = LinearSolver::solve(K, -R);
    VectorX du_bar2 = LinearSolver::solve(K, -q);

    Real denom_corr = v0.dot(du_bar2) + options.psi * options.psi * delta_lambda_predictor;
    if (std::abs(denom_corr) < 1e-13) {
      throw std::runtime_error(
          "ArcLengthSolver: degenerate correction (normal-plane direction nearly parallel to "
          "the residual's lambda-sensitivity)");
    }
    Real d_delta_lambda = -(v0.dot(du_bar1)) / denom_corr;
    VectorX d_delta_u = du_bar1 + d_delta_lambda * du_bar2;

    delta_u += d_delta_u;
    delta_lambda += d_delta_lambda;
  }

  result.delta_u = delta_u;
  result.delta_lambda = delta_lambda;
  result.iterations = iter;
  return result;
}

}  // namespace

ArcLengthSolver::ArcLengthSolver() : ArcLengthSolver(Options{}) {}

ArcLengthSolver::ArcLengthSolver(Options options) : options_(options) {}

std::vector<ArcLengthSolver::StepResult> ArcLengthSolver::solve(const ResidualFunction& residual,
                                                                 const TangentFunction& tangent,
                                                                 const VectorX& u0,
                                                                 Real lambda0) const {
  std::vector<StepResult> path;

  VectorX u_current = u0;
  Real lambda_current = lambda0;

  // "Previous step's increment", used for the sign-continuity
  // criterion. There is no real previous step yet, so this just
  // encodes the requested starting direction via its lambda component.
  VectorX prev_delta_u = VectorX::Zero(u0.size());
  Real prev_delta_lambda = options_.initial_direction;

  for (int step = 0; step < options_.num_steps; ++step) {
    Real arc_length = options_.arc_length;
    AttemptResult attempt;

    for (int retry = 0; retry <= options_.max_retries_per_step; ++retry) {
      attempt = attempt_step(residual, tangent, options_, u_current, lambda_current, prev_delta_u,
                              prev_delta_lambda, arc_length);
      if (attempt.converged) {
        break;
      }
      arc_length *= options_.retry_shrink_factor;  // smaller step, try again
    }

    StepResult result;
    result.u = u_current + attempt.delta_u;
    result.lambda = lambda_current + attempt.delta_lambda;
    result.converged = attempt.converged;
    result.iterations = attempt.iterations;
    path.push_back(result);

    if (!attempt.converged) {
      break;  // every retry failed -- report what converged, same policy as elsewhere
    }

    prev_delta_u = attempt.delta_u;
    prev_delta_lambda = attempt.delta_lambda;
    u_current = result.u;
    lambda_current = result.lambda;
  }

  return path;
}

}  // namespace foam
