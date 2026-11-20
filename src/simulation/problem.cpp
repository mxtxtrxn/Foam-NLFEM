#include "foam/simulation/problem.hpp"

#include <stdexcept>
#include <utility>

namespace foam {

Problem::Problem(Mesh mesh, const MaterialModel& material,
                  const FiniteElement& element_formulation, DirichletBCs dirichlet_bcs,
                  NeumannBCs neumann_bcs)
    : mesh_(std::move(mesh)),
      material_(material),
      element_formulation_(element_formulation),
      dirichlet_bcs_(std::move(dirichlet_bcs)),
      neumann_bcs_(std::move(neumann_bcs)) {}

std::vector<std::pair<Real, Real>> Problem::run_compression_test(int /*num_steps*/) {
  // TODO(roadmap step 2): for each load step,
  //   1. scale dirichlet_bcs_/neumann_bcs_ targets by step fraction
  //   2. build residual/tangent lambdas around an Assembler
  //      (capturing the current step's BC targets)
  //   3. call NewtonRaphsonSolver::solve(...)
  //   4. record (nominal strain, nominal stress) from the reaction force
  //      at the constrained face and the current specimen height
  throw std::logic_error("Problem::run_compression_test: not implemented yet");
}

}  // namespace foam
