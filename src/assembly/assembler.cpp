#include "foam/assembly/assembler.hpp"

namespace foam {

Assembler::Assembler(const Mesh& mesh, const FiniteElement& element_formulation,
                      const MaterialModel& material, int dofs_per_node)
    : mesh_(mesh),
      element_formulation_(element_formulation),
      material_(material),
      dofs_per_node_(dofs_per_node) {}

namespace {

// Gathers one element's reference positions and current nodal
// displacements out of the global arrays, as num_nodes x dofs_per_node
// matrices (the shape every FiniteElement method expects). Relies on
// Mesh::nodes()[id] being indexed by that same id -- true by
// construction, since Mesh::add_node hands out ids as nodes_.size() in
// call order.
void gather_element_data(const Mesh& mesh, const ElementConnectivity& elem, int dofs_per_node,
                          const VectorX& global_displacements, MatrixX& reference_positions,
                          MatrixX& nodal_displacements) {
  int num_nodes = static_cast<int>(elem.node_ids.size());
  reference_positions.resize(num_nodes, dofs_per_node);
  nodal_displacements.resize(num_nodes, dofs_per_node);
  for (int a = 0; a < num_nodes; ++a) {
    std::size_t node_id = elem.node_ids[a];
    const Vector3& X = mesh.nodes()[node_id].reference_position();
    for (int d = 0; d < dofs_per_node; ++d) {
      reference_positions(a, d) = X(d);
      nodal_displacements(a, d) = global_displacements(node_id * dofs_per_node + d);
    }
  }
}

}  // namespace

VectorX Assembler::assemble_residual(const VectorX& global_displacements) const {
  VectorX result = VectorX::Zero(mesh_.num_dofs(dofs_per_node_));

  MatrixX reference_positions;
  MatrixX nodal_displacements;
  for (const auto& elem : mesh_.elements()) {
    gather_element_data(mesh_, elem, dofs_per_node_, global_displacements, reference_positions,
                         nodal_displacements);

    VectorX f_e =
        element_formulation_.internal_force(reference_positions, nodal_displacements, material_);

    int num_nodes = static_cast<int>(elem.node_ids.size());
    for (int a = 0; a < num_nodes; ++a) {
      std::size_t node_id = elem.node_ids[a];
      for (int d = 0; d < dofs_per_node_; ++d) {
        result(node_id * dofs_per_node_ + d) += f_e(a * dofs_per_node_ + d);
      }
    }
  }
  return result;
}

SparseMatrix Assembler::assemble_tangent(const VectorX& global_displacements) const {
  std::size_t total_dofs = mesh_.num_dofs(dofs_per_node_);
  std::vector<Triplet> triplets;

  MatrixX reference_positions;
  MatrixX nodal_displacements;
  for (const auto& elem : mesh_.elements()) {
    gather_element_data(mesh_, elem, dofs_per_node_, global_displacements, reference_positions,
                         nodal_displacements);

    MatrixX K_e = element_formulation_.tangent_stiffness(reference_positions,
                                                          nodal_displacements, material_);

    int num_nodes = static_cast<int>(elem.node_ids.size());
    for (int a = 0; a < num_nodes; ++a) {
      std::size_t node_a = elem.node_ids[a];
      for (int da = 0; da < dofs_per_node_; ++da) {
        std::size_t global_row = node_a * dofs_per_node_ + da;
        for (int b = 0; b < num_nodes; ++b) {
          std::size_t node_b = elem.node_ids[b];
          for (int db = 0; db < dofs_per_node_; ++db) {
            std::size_t global_col = node_b * dofs_per_node_ + db;
            Real value = K_e(a * dofs_per_node_ + da, b * dofs_per_node_ + db);
            triplets.emplace_back(static_cast<int>(global_row), static_cast<int>(global_col),
                                   value);
          }
        }
      }
    }
  }

  SparseMatrix K(static_cast<int>(total_dofs), static_cast<int>(total_dofs));
  K.setFromTriplets(triplets.begin(), triplets.end());
  return K;
}

}  // namespace foam
