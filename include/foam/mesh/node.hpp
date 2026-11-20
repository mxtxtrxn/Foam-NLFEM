#pragma once

#include <cstddef>

#include "foam/core/types.hpp"

namespace foam {

// A single mesh node: an id plus its reference (undeformed) coordinates.
// Current (deformed) coordinates are NOT stored here -- they live in the
// global displacement vector owned by Problem/Assembler, so Node stays
// a lightweight, immutable geometric fact about the mesh.
class Node {
 public:
  Node(std::size_t id, Vector3 reference_position);

  std::size_t id() const { return id_; }
  const Vector3& reference_position() const { return reference_position_; }

 private:
  std::size_t id_;
  Vector3 reference_position_;
};

}  // namespace foam
