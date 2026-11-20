#include "foam/assembly/assembler.hpp"

#include <stdexcept>

namespace foam {

Assembler::Assembler(const Mesh& mesh, const FiniteElement& element_formulation,
                      const MaterialModel& material, int dofs_per_node)
    : mesh_(mesh),
      element_formulation_(element_formulation),
      material_(material),
      dofs_per_node_(dofs_per_node) {}

VectorX Assembler::assemble_residual(const VectorX& /*global_displacements*/) const {
  // TODO(roadmap step 2): for each element in mesh_, gather its nodal
  // reference positions + displacements, call
  // element_formulation_.internal_force(...), and scatter-add into the
  // global residual using each node's global DOF indices
  // (node_id * dofs_per_node_ + local_dof).
  throw std::logic_error("Assembler::assemble_residual: not implemented yet");
}

SparseMatrix Assembler::assemble_tangent(const VectorX& /*global_displacements*/) const {
  // TODO(roadmap step 2): same loop as assemble_residual, but building a
  // list of Eigen::Triplet entries from tangent_stiffness() and calling
  // SparseMatrix::setFromTriplets at the end.
  throw std::logic_error("Assembler::assemble_tangent: not implemented yet");
}

}  // namespace foam
