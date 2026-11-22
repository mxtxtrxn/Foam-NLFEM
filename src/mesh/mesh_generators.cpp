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

Mesh make_structured_grid_2d(Real width, Real height, int nx, int ny) {
  if (nx < 1 || ny < 1) {
    throw std::invalid_argument("make_structured_grid_2d: nx and ny must be >= 1");
  }

  Mesh mesh;

  // Node (i,j), i=0..nx, j=0..ny, at (i*dx, j*dy, 0), indexed row-major
  // as j*(nx+1)+i. This indexing is relied on below when building each
  // element's connectivity, and gives the id Mesh::add_node assigns
  // (ids are handed out in call order, 0, 1, 2, ...).
  Real dx = width / nx;
  Real dy = height / ny;
  for (int j = 0; j <= ny; ++j) {
    for (int i = 0; i <= nx; ++i) {
      mesh.add_node(Vector3(i * dx, j * dy, 0.0));
    }
  }

  auto node_id = [nx](int i, int j) { return static_cast<std::size_t>(j * (nx + 1) + i); };

  // Bottom-left, bottom-right, top-right, top-left (CCW), matching the
  // Quad4 shape-function node ordering Quad4Element uses.
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      mesh.add_element(ElementType::kQuad4, {node_id(i, j), node_id(i + 1, j),
                                              node_id(i + 1, j + 1), node_id(i, j + 1)});
    }
  }

  return mesh;
}

}  // namespace foam
