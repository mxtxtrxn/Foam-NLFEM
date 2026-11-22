#pragma once

#include <functional>
#include <vector>

#include "foam/core/types.hpp"

namespace foam {

// Arc-length (Riks/Wempner "normal plane" variant) continuation. Plain
// Newton-Raphson under load OR displacement control gets stuck at
// limit points -- any point where the tangent K = dR/du becomes
// singular, which is exactly what a snap-through/softening instability
// looks like. Arc-length traces the equilibrium path by treating the
// control parameter lambda as an extra unknown, solved for alongside u
// at each step, subject to a constraint on the combined step length.
//
// DESIGN NOTE (revised from the Step-1-era declaration): this version
// takes R(u, lambda) as a single abstract function and computes
// dR/dlambda internally via finite differences, rather than requiring
// a separate "reference load vector" assuming R(u,lambda) = F_int(u) -
// lambda*F_ref (pure load control). That assumption doesn't fit
// Problem::run_compression_test()'s displacement-control setup, where
// lambda is the Dirichlet-BC ramp fraction t and dR/dlambda depends on
// the free-constrained coupling of the tangent stiffness, not a fixed
// vector. Computing dR/dlambda generically, the same way
// numerical_tangent_moduli() computes a material tangent, works for
// either case without the solver needing to know which one it's in.
class ArcLengthSolver {
 public:
  struct Options {
    Real arc_length = 0.1;           // target combined step length (sqrt(||du||^2 + psi^2*dlambda^2))
    Real psi = 1.0;                   // weight of dlambda^2 in the arc-length constraint
    int max_iterations_per_step = 30;
    int num_steps = 100;
    Real tolerance = 1e-8;             // ABSOLUTE residual-norm tolerance (as in NewtonRaphsonSolver)
    Real lambda_perturbation = 1e-6;   // step for the internal finite-difference dR/dlambda
    Real initial_direction = 1.0;      // sign only (+1 or -1): which way to start lambda moving
    int max_retries_per_step = 4;      // halve the arc length and retry this many times on failure
    Real retry_shrink_factor = 0.5;
  };

  struct StepResult {
    VectorX u;
    Real lambda = 0.0;
    bool converged = false;
    int iterations = 0;
  };

  using ResidualFunction = std::function<VectorX(const VectorX& u, Real lambda)>;
  using TangentFunction = std::function<SparseMatrix(const VectorX& u, Real lambda)>;

  ArcLengthSolver();
  explicit ArcLengthSolver(Options options);

  // Traces the equilibrium path starting from the (presumed converged)
  // point (u0, lambda0), taking up to options.num_steps arc-length
  // steps. Returns one StepResult per step attempted; stops (returning
  // fewer than num_steps results) at the first step whose corrector
  // iterations fail to converge, rather than continuing from an
  // unreliable point -- the same policy used elsewhere in this project
  // (foam::solve_uniaxial_point, Problem::run_compression_test).
  std::vector<StepResult> solve(const ResidualFunction& residual, const TangentFunction& tangent,
                                 const VectorX& u0, Real lambda0) const;

 private:
  Options options_;
};

}  // namespace foam
