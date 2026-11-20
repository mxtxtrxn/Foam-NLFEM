#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace foam {

// Pure mesh-level connectivity: which node ids make up this element and
// what kind of element it is. This is NOT the finite-element formulation
// (shape functions, quadrature, stiffness) -- that lives in foam/elements.
// Keeping topology and formulation separate means the mesh module has no
// numerical-methods knowledge, and the element module has no I/O/file
// -format knowledge.
enum class ElementType { kQuad4, kHex8, kAxisymmetricQuad4 };

inline std::string to_string(ElementType type) {
  switch (type) {
    case ElementType::kQuad4:
      return "Quad4";
    case ElementType::kHex8:
      return "Hex8";
    case ElementType::kAxisymmetricQuad4:
      return "AxisymmetricQuad4";
  }
  return "Unknown";
}

struct ElementConnectivity {
  std::size_t id;
  ElementType type;
  std::vector<std::size_t> node_ids;
};

}  // namespace foam
