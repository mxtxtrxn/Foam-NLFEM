#pragma once

#include <vector>

#include "foam/mesh/element.hpp"
#include "foam/mesh/node.hpp"

namespace foam {

// Owns the geometric description of a problem: nodes + element
// connectivity. Mesh knows nothing about materials, boundary conditions,
// or solvers -- those are composed on top of it in foam::Problem.
class Mesh {
 public:
  std::size_t add_node(Vector3 reference_position);
  std::size_t add_element(ElementType type, std::vector<std::size_t> node_ids);

  const std::vector<Node>& nodes() const { return nodes_; }
  const std::vector<ElementConnectivity>& elements() const { return elements_; }

  std::size_t num_nodes() const { return nodes_.size(); }
  std::size_t num_elements() const { return elements_.size(); }

  // Each node has `dofs_per_node` degrees of freedom (e.g. 2 for 2D
  // displacement, 3 for 3D). Total system size = num_nodes * dofs_per_node.
  std::size_t num_dofs(std::size_t dofs_per_node) const {
    return nodes_.size() * dofs_per_node;
  }

 private:
  std::vector<Node> nodes_;
  std::vector<ElementConnectivity> elements_;
};

}  // namespace foam
