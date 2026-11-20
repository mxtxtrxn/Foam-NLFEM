#include "foam/solvers/arc_length.hpp"

#include <stdexcept>

namespace foam {

ArcLengthSolver::ArcLengthSolver(Options options) : options_(options) {}

void ArcLengthSolver::solve(const ResidualFunction& /*residual*/,
                             const TangentFunction& /*tangent*/,
                             const VectorX& /*reference_load_vector*/) const {
  // TODO(roadmap step 3): spherical or cylindrical arc-length constraint,
  // predictor-corrector with the load factor as an extra unknown.
  throw std::logic_error("ArcLengthSolver::solve: not implemented yet");
}

}  // namespace foam
