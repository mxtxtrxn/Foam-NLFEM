#pragma once

#include "foam/mesh/mesh.hpp"

namespace foam {

// Builds a single unit Hex8 element (a 1x1x1 cube). This is the geometry
// for the MVP step: validate a material model (e.g. Ogden hyperfoam)
// under a homogeneous deformation before any real assembly is involved.
Mesh make_single_hex_element();

// Builds a structured 2D grid of Quad4 elements over [0, width] x [0, height].
// nx, ny = number of elements along each direction.
// Deferred to the "2D/axisymmetric specimen" stage of the roadmap.
Mesh make_structured_grid_2d(Real width, Real height, int nx, int ny);

}  // namespace foam
