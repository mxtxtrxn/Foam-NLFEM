#pragma once

#include "foam/core/types.hpp"

namespace foam {

// Thin wrapper over Eigen's sparse direct solver. Kept as its own small
// class (rather than calling Eigen::SparseLU directly from
// NewtonRaphsonSolver) so the linear-algebra backend can be swapped later
// (e.g. an iterative solver + preconditioner for very large 3D meshes)
// without touching the nonlinear solvers built on top of it.
class LinearSolver {
 public:
  // Solves K * x = rhs for x. Throws std::runtime_error if K is singular
  // or factorization otherwise fails.
  static VectorX solve(const SparseMatrix& K, const VectorX& rhs);
};

}  // namespace foam
