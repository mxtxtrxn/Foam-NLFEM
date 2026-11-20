#pragma once

#include <cstddef>
#include <vector>

#include "foam/core/types.hpp"

namespace foam {

// A single prescribed-displacement constraint: node `node_id`, local
// degree of freedom `dof` (0=x, 1=y, 2=z), prescribed to `value`.
// These are plain data; applying them (penalty, elimination, or
// Lagrange-multiplier method) is the solver's job, not this struct's.
struct DirichletBC {
  std::size_t node_id;
  int dof;
  Real value;
};

using DirichletBCs = std::vector<DirichletBC>;

}  // namespace foam
