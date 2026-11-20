#include "foam/solvers/linear_solver.hpp"

#include <Eigen/SparseLU>
#include <stdexcept>

namespace foam {

VectorX LinearSolver::solve(const SparseMatrix& K, const VectorX& rhs) {
  Eigen::SparseLU<SparseMatrix> factorization;
  factorization.analyzePattern(K);
  factorization.factorize(K);
  if (factorization.info() != Eigen::Success) {
    throw std::runtime_error("LinearSolver::solve: factorization failed (singular K?)");
  }
  VectorX x = factorization.solve(rhs);
  if (factorization.info() != Eigen::Success) {
    throw std::runtime_error("LinearSolver::solve: back-substitution failed");
  }
  return x;
}

}  // namespace foam
