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

TEST_CASE("make_structured_grid_2d: correct node/element counts and corner positions") {
  auto mesh = foam::make_structured_grid_2d(/*width=*/2.0, /*height=*/1.0, /*nx=*/2, /*ny=*/3);

  CHECK(mesh.num_nodes() == 3 * 4);  // (nx+1)*(ny+1)
  CHECK(mesh.num_elements() == 2 * 3);  // nx*ny

  // Corner nodes, by construction node_id(i,j) = j*(nx+1)+i.
  CHECK(mesh.nodes()[0].reference_position().isApprox(foam::Vector3(0.0, 0.0, 0.0)));  // (0,0)
  CHECK(mesh.nodes()[2].reference_position().isApprox(foam::Vector3(2.0, 0.0, 0.0)));  // (nx,0)
  CHECK(mesh.nodes()[9].reference_position().isApprox(foam::Vector3(0.0, 1.0, 0.0)));  // (0,ny)
  CHECK(mesh.nodes()[11].reference_position().isApprox(foam::Vector3(2.0, 1.0, 0.0)));  // (nx,ny)

  // All elements are Quad4 with 4 distinct node ids, CCW starting
  // bottom-left (matches Quad4Element's shape-function convention).
  for (const auto& elem : mesh.elements()) {
    CHECK(elem.type == foam::ElementType::kQuad4);
    REQUIRE(elem.node_ids.size() == 4);
    for (int a = 0; a < 4; ++a) {
      for (int b = a + 1; b < 4; ++b) {
        CHECK(elem.node_ids[a] != elem.node_ids[b]);
      }
    }
  }
}

TEST_CASE("make_structured_grid_2d: rejects nx or ny < 1") {
  CHECK_THROWS_AS(foam::make_structured_grid_2d(1.0, 1.0, 0, 2), std::invalid_argument);
  CHECK_THROWS_AS(foam::make_structured_grid_2d(1.0, 1.0, 2, 0), std::invalid_argument);
}
