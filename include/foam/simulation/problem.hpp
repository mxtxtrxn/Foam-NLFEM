#pragma once

#include <vector>

#include "foam/boundary_conditions/dirichlet_bc.hpp"
#include "foam/boundary_conditions/neumann_bc.hpp"
#include "foam/elements/element_base.hpp"
#include "foam/materials/material_base.hpp"
#include "foam/mesh/mesh.hpp"

namespace foam {

// The highest-level object: wires a Mesh + MaterialModel + FiniteElement
// formulation + boundary conditions into something that can be driven
// through a load/displacement history and produce a result (e.g. a
// stress-strain curve). This is deliberately the LAST thing to implement
// -- everything it depends on (Assembler, NewtonRaphsonSolver, etc.)
// needs to work and be unit-tested in isolation first.
class Problem {
 public:
  Problem(Mesh mesh, const MaterialModel& material, const FiniteElement& element_formulation,
          DirichletBCs dirichlet_bcs, NeumannBCs neumann_bcs);

  // Applies `num_steps` increments of the prescribed Dirichlet/Neumann
  // BCs (scaled linearly from 0 to 1), solving the nonlinear system at
  // each step, and returns the (nominal strain, nominal stress) pairs
  // useful for a compression stress-strain curve.
  //
  // TODO(roadmap step 2, once Assembler + NewtonRaphsonSolver are wired
  // together here): implement the load-stepping loop.
  std::vector<std::pair<Real, Real>> run_compression_test(int num_steps);

 private:
  Mesh mesh_;
  const MaterialModel& material_;
  const FiniteElement& element_formulation_;
  DirichletBCs dirichlet_bcs_;
  NeumannBCs neumann_bcs_;
};

}  // namespace foam
