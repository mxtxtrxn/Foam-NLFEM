#pragma once

#include <cstddef>
#include <vector>

#include "foam/core/types.hpp"

namespace foam {

// A concentrated nodal force. Distributed (surface pressure) loads can be
// added later as a separate type once an element has face/edge
// integration -- keeping this struct minimal for now.
struct NeumannBC {
  std::size_t node_id;
  int dof;
  Real value;
};

using NeumannBCs = std::vector<NeumannBC>;

}  // namespace foam
