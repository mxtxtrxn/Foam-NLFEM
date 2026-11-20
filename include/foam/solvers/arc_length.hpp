#pragma once

#include <functional>

#include "foam/core/types.hpp"

namespace foam {

// Arc-length (Riks-type) continuation. Plain Newton-Raphson under
// load/displacement control gets stuck at limit points -- exactly what
// happens at the onset of the foam's collapse plateau, where the
// load-displacement curve can go flat or even turn over locally. This
// solver follows the equilibrium path by constraining a combination of
// the load factor and displacement increment, not just displacement.
//
// Deferred to roadmap step 3: get plain Newton-Raphson working end to
// end on the 2D specimen first; add this once the plateau/instability
// region is actually being reached and plain NR starts failing there.
class ArcLengthSolver {
 public:
  struct Options {
    Real initial_arc_length = 0.1;
    int max_iterations_per_step = 30;
    int num_steps = 100;
    Real tolerance = 1e-8;
  };

  using ResidualFunction = std::function<VectorX(const VectorX& displacement, Real load_factor)>;
  using TangentFunction = std::function<SparseMatrix(const VectorX& displacement, Real load_factor)>;

  explicit ArcLengthSolver(Options options);

  // TODO(roadmap step 3): implement the Riks predictor-corrector loop.
  void solve(const ResidualFunction& residual, const TangentFunction& tangent,
             const VectorX& reference_load_vector) const;

 private:
  Options options_;
};

}  // namespace foam
