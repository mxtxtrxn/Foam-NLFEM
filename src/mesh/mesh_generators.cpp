#include "foam/mesh/mesh_generators.hpp"

#include <stdexcept>

namespace foam {

Mesh make_single_hex_element() {
  Mesh mesh;
  // Standard Hex8 node ordering: bottom face CCW, then top face CCW.
  std::size_t n0 = mesh.add_node(Vector3(0, 0, 0));
  std::size_t n1 = mesh.add_node(Vector3(1, 0, 0));
  std::size_t n2 = mesh.add_node(Vector3(1, 1, 0));
  std::size_t n3 = mesh.add_node(Vector3(0, 1, 0));
  std::size_t n4 = mesh.add_node(Vector3(0, 0, 1));
  std::size_t n5 = mesh.add_node(Vector3(1, 0, 1));
  std::size_t n6 = mesh.add_node(Vector3(1, 1, 1));
  std::size_t n7 = mesh.add_node(Vector3(0, 1, 1));

  mesh.add_element(ElementType::kHex8, {n0, n1, n2, n3, n4, n5, n6, n7});
  return mesh;
}

Mesh make_structured_grid_2d(Real /*width*/, Real /*height*/, int /*nx*/, int /*ny*/) {
  // TODO(roadmap step 2): generate a structured Quad4 grid for the
  // 2D / axisymmetric compression-specimen case.
  throw std::logic_error("make_structured_grid_2d: not implemented yet");
}

}  // namespace foam
