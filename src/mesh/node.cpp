#include "foam/mesh/node.hpp"

namespace foam {

Node::Node(std::size_t id, Vector3 reference_position)
    : id_(id), reference_position_(std::move(reference_position)) {}

}  // namespace foam
