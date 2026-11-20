#include "doctest.h"
#include "foam/mesh/mesh.hpp"
#include "foam/mesh/mesh_generators.hpp"

TEST_CASE("Mesh: add_node/add_element basic bookkeeping") {
  foam::Mesh mesh;
  auto n0 = mesh.add_node(foam::Vector3(0, 0, 0));
  auto n1 = mesh.add_node(foam::Vector3(1, 0, 0));
  CHECK(mesh.num_nodes() == 2);
  CHECK(n0 == 0);
  CHECK(n1 == 1);

  auto e0 = mesh.add_element(foam::ElementType::kQuad4, {n0, n1, n0, n1});
  CHECK(mesh.num_elements() == 1);
  CHECK(e0 == 0);
  CHECK(mesh.num_dofs(/*dofs_per_node=*/2) == 4);
}

TEST_CASE("make_single_hex_element produces 8 nodes and one Hex8 element") {
  auto mesh = foam::make_single_hex_element();
  CHECK(mesh.num_nodes() == 8);
  CHECK(mesh.num_elements() == 1);
  CHECK(mesh.elements()[0].type == foam::ElementType::kHex8);
  CHECK(mesh.elements()[0].node_ids.size() == 8);

  // Sanity check on geometry: it should be the unit cube [0,1]^3.
  for (const auto& node : mesh.nodes()) {
    const auto& p = node.reference_position();
    CHECK(p.x() >= 0.0);
    CHECK(p.x() <= 1.0);
    CHECK(p.y() >= 0.0);
    CHECK(p.y() <= 1.0);
    CHECK(p.z() >= 0.0);
    CHECK(p.z() <= 1.0);
  }
}

TEST_CASE("make_structured_grid_2d is a documented not-yet-implemented stub") {
  CHECK_THROWS_AS(foam::make_structured_grid_2d(1.0, 1.0, 2, 2), std::logic_error);
}
