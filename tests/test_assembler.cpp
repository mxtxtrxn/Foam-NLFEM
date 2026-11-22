#include "doctest.h"
#include "foam/assembly/assembler.hpp"
#include "foam/elements/quad4.hpp"
#include "foam/materials/neo_hookean.hpp"
#include "foam/mesh/mesh_generators.hpp"

TEST_CASE("Assembler: zero displacement everywhere gives zero global internal force") {
  auto mesh = foam::make_structured_grid_2d(2.0, 2.0, 3, 3);
  foam::NeoHookean material(1.0e5, 1.5e5);
  foam::Quad4Element element;
  foam::Assembler assembler(mesh, element, material, /*dofs_per_node=*/2);

  foam::VectorX disp = foam::VectorX::Zero(mesh.num_dofs(2));
  foam::VectorX f = assembler.assemble_residual(disp);
  CHECK(f.norm() < 1e-6);
}

TEST_CASE("Assembler: PATCH TEST -- a homogeneous deformation applied everywhere gives zero "
          "internal force at every INTERIOR node") {
  // The classic FEM patch test: a uniform stress/strain field satisfies
  // equilibrium identically (no body force), so at any node not on the
  // mesh boundary -- where contributions from every surrounding element
  // must cancel exactly -- the assembled internal force must vanish.
  // This is a strong, implementation-independent check of the
  // element-to-global scatter/assembly logic, on a mesh with actual
  // interior nodes (unlike the single-element tests in test_quad4.cpp).
  int nx = 3, ny = 3;
  auto mesh = foam::make_structured_grid_2d(2.0, 2.0, nx, ny);
  foam::NeoHookean material(1.0e5, 1.5e5);
  foam::Quad4Element element;
  foam::Assembler assembler(mesh, element, material, /*dofs_per_node=*/2);

  foam::Real axial = 0.85;
  foam::Real lateral = 1.05;
  foam::VectorX disp = foam::VectorX::Zero(mesh.num_dofs(2));
  for (const auto& node : mesh.nodes()) {
    const foam::Vector3& X = node.reference_position();
    disp(node.id() * 2 + 0) = (axial - 1.0) * X(0);
    disp(node.id() * 2 + 1) = (lateral - 1.0) * X(1);
  }

  foam::VectorX f = assembler.assemble_residual(disp);

  // Interior nodes: i in [1, nx-1], j in [1, ny-1] (strictly inside,
  // surrounded by elements on all four sides).
  int interior_count = 0;
  for (int j = 1; j < ny; ++j) {
    for (int i = 1; i < nx; ++i) {
      std::size_t node_id = static_cast<std::size_t>(j * (nx + 1) + i);
      CHECK(std::abs(f(node_id * 2 + 0)) < 1e-6);
      CHECK(std::abs(f(node_id * 2 + 1)) < 1e-6);
      ++interior_count;
    }
  }
  REQUIRE(interior_count == (nx - 1) * (ny - 1));
}

TEST_CASE("Assembler: assemble_tangent produces a symmetric global stiffness matrix") {
  // Each element's tangent_stiffness is symmetric (material part is
  // B_L^T D B_L with symmetric D; geometric part is symmetric by
  // construction), and assembly is a linear scatter, so the global
  // matrix must be symmetric too.
  auto mesh = foam::make_structured_grid_2d(1.5, 1.0, 2, 2);
  foam::NeoHookean material(1.0e5, 1.5e5);
  foam::Quad4Element element;
  foam::Assembler assembler(mesh, element, material, 2);

  foam::VectorX disp = foam::VectorX::Zero(mesh.num_dofs(2));
  for (const auto& node : mesh.nodes()) {
    const foam::Vector3& X = node.reference_position();
    disp(node.id() * 2 + 0) = 0.05 * X(0) - 0.02 * X(1);
    disp(node.id() * 2 + 1) = -0.03 * X(1);
  }

  foam::SparseMatrix K = assembler.assemble_tangent(disp);
  foam::MatrixX K_dense = foam::MatrixX(K);
  CHECK((K_dense - K_dense.transpose()).cwiseAbs().maxCoeff() < 1e-6);
}
