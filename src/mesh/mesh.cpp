#include "foam/mesh/mesh.hpp"

namespace foam {

std::size_t Mesh::add_node(Vector3 reference_position) {
  std::size_t id = nodes_.size();
  nodes_.emplace_back(id, std::move(reference_position));
  return id;
}

std::size_t Mesh::add_element(ElementType type, std::vector<std::size_t> node_ids) {
  std::size_t id = elements_.size();
  elements_.push_back(ElementConnectivity{id, type, std::move(node_ids)});
  return id;
}

}  // namespace foam
