#pragma once

#include "foam/elements/element_base.hpp"
#include "foam/mesh/mesh.hpp"

namespace foam {

// Loops over all elements in a Mesh, calls each FiniteElement's
// internal_force()/tangent_stiffness(), and scatters the element-level
// contributions into the global residual vector / global sparse tangent
// matrix. This is the one place that knows about global DOF numbering;
// everything below it (elements, materials) only ever sees local,
// per-element data.
class Assembler {
 public:
  Assembler(const Mesh& mesh, const FiniteElement& element_formulation,
            const MaterialModel& material, int dofs_per_node);

  // Global residual: R(u) = F_internal(u) - F_external. Returned vector
  // has size mesh.num_dofs(dofs_per_node).
  VectorX assemble_residual(const VectorX& global_displacements) const;

  // Global tangent stiffness, assembled as a sparse matrix (most entries
  // are zero since each node only couples to its neighboring elements).
  SparseMatrix assemble_tangent(const VectorX& global_displacements) const;

 private:
  const Mesh& mesh_;
  const FiniteElement& element_formulation_;
  const MaterialModel& material_;
  int dofs_per_node_;
};

}  // namespace foam
