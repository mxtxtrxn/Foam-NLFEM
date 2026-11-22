#include "foam/simulation/problem.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "foam/assembly/assembler.hpp"
#include "foam/solvers/newton_raphson.hpp"

namespace foam {

Problem::Problem(Mesh mesh, const MaterialModel& material,
                  const FiniteElement& element_formulation, DirichletBCs dirichlet_bcs,
                  NeumannBCs neumann_bcs)
    : mesh_(std::move(mesh)),
      material_(material),
      element_formulation_(element_formulation),
      dirichlet_bcs_(std::move(dirichlet_bcs)),
      neumann_bcs_(std::move(neumann_bcs)) {}

std::vector<std::pair<Real, Real>> Problem::run_compression_test(int num_steps) {
  // This method is specific to the use case its name promises (a
  // uniaxial compression test on a rectangular specimen), with these
  // conventions, rather than being a fully general BC-driven solver:
  //   - Y is the loading (axial) direction; X is lateral.
  //   - dirichlet_bcs_ entries on dof=1 (y) with a nonzero target are
  //     "the loaded edge" (assumed a single flat rigid platen: all such
  //     targets must be equal); their value is the FULL-load
  //     displacement, ramped linearly from 0 over num_steps.
  //   - Reference height/width come from the mesh's own bounding box
  //     (max Y, max X among nodes), i.e. the specimen is assumed to
  //     occupy [0,width] x [0,height] in the reference configuration.
  //   - Reported "stress" is the summed reaction force at the loaded
  //     edge divided by the reference width (unit thickness, plane
  //     strain) -- a nominal (1st Piola-Kirchhoff-like) stress, the
  //     plane-strain analogue of UniaxialPoint::nominal_stress_axial in
  //     foam/material_point/uniaxial_test.hpp.
  if (neumann_bcs_.empty() == false) {
    throw std::logic_error(
        "Problem::run_compression_test: Neumann BCs are not supported by this method yet "
        "(it is a purely displacement-controlled test)");
  }

  int dofs_per_node = element_formulation_.num_dofs_per_node();
  if (dofs_per_node != 2) {
    throw std::logic_error(
        "Problem::run_compression_test: only 2D (dofs_per_node == 2) elements are supported");
  }
  std::size_t total_dofs = mesh_.num_dofs(dofs_per_node);

  Assembler assembler(mesh_, element_formulation_, material_, dofs_per_node);

  // Constrained-dof bookkeeping: target value is the FULL-load (step
  // fraction t=1) displacement; applied value at step t is t*target.
  std::vector<char> is_constrained(total_dofs, 0);
  std::vector<Real> target_value(total_dofs, 0.0);
  for (const auto& bc : dirichlet_bcs_) {
    std::size_t gdof = bc.node_id * dofs_per_node + bc.dof;
    if (gdof >= total_dofs) {
      throw std::out_of_range("Problem::run_compression_test: a DirichletBC references a dof "
                               "outside the mesh");
    }
    is_constrained[gdof] = 1;
    target_value[gdof] = bc.value;
  }

  std::vector<std::size_t> free_dofs;
  std::vector<int> global_to_free(total_dofs, -1);
  for (std::size_t g = 0; g < total_dofs; ++g) {
    if (!is_constrained[g]) {
      global_to_free[g] = static_cast<int>(free_dofs.size());
      free_dofs.push_back(g);
    }
  }
  std::size_t num_free = free_dofs.size();
  if (num_free == 0) {
    throw std::logic_error("Problem::run_compression_test: every dof is constrained (nothing "
                            "to solve)");
  }

  // "Loaded" dofs: y-dof Dirichlet entries with a nonzero target.
  std::vector<std::size_t> loaded_dofs;
  Real applied_target = 0.0;
  for (const auto& bc : dirichlet_bcs_) {
    if (bc.dof == 1 && std::abs(bc.value) > 1e-12) {
      std::size_t gdof = bc.node_id * dofs_per_node + bc.dof;
      loaded_dofs.push_back(gdof);
      if (loaded_dofs.size() == 1) {
        applied_target = bc.value;
      } else if (std::abs(bc.value - applied_target) > 1e-9 * std::abs(applied_target)) {
        throw std::logic_error(
            "Problem::run_compression_test: loaded-edge y-displacement targets are not all "
            "equal (expected a single flat rigid platen)");
      }
    }
  }
  if (loaded_dofs.empty()) {
    throw std::logic_error(
        "Problem::run_compression_test: no nonzero y-displacement Dirichlet BC found -- "
        "nothing is being loaded");
  }

  Real height = 0.0;
  Real width = 0.0;
  for (const auto& node : mesh_.nodes()) {
    height = std::max(height, node.reference_position()(1));
    width = std::max(width, node.reference_position()(0));
  }
  if (height <= 0.0 || width <= 0.0) {
    throw std::logic_error("Problem::run_compression_test: degenerate mesh bounding box");
  }

  auto build_full_displacement = [&](Real t, const VectorX& free_values) {
    VectorX full = VectorX::Zero(static_cast<Eigen::Index>(total_dofs));
    for (std::size_t g = 0; g < total_dofs; ++g) {
      full(static_cast<Eigen::Index>(g)) =
          is_constrained[g] ? t * target_value[g] : free_values(global_to_free[g]);
    }
    return full;
  };

  std::vector<std::pair<Real, Real>> curve;
  curve.reserve(static_cast<std::size_t>(num_steps) + 1);
  curve.emplace_back(0.0, 0.0);  // reference (undeformed) state

  VectorX free_guess = VectorX::Zero(static_cast<Eigen::Index>(num_free));

  for (int step = 1; step <= num_steps; ++step) {
    Real t = static_cast<Real>(step) / static_cast<Real>(num_steps);

    auto residual = [&](const VectorX& free_values) -> VectorX {
      VectorX full = build_full_displacement(t, free_values);
      VectorX f_int_full = assembler.assemble_residual(full);
      VectorX r(static_cast<Eigen::Index>(num_free));
      for (std::size_t k = 0; k < num_free; ++k) {
        r(static_cast<Eigen::Index>(k)) = f_int_full(static_cast<Eigen::Index>(free_dofs[k]));
      }
      return r;
    };

    auto tangent = [&](const VectorX& free_values) -> SparseMatrix {
      VectorX full = build_full_displacement(t, free_values);
      SparseMatrix K_full = assembler.assemble_tangent(full);
      std::vector<Triplet> triplets;
      for (int k = 0; k < K_full.outerSize(); ++k) {
        for (SparseMatrix::InnerIterator it(K_full, k); it; ++it) {
          int gr = global_to_free[static_cast<std::size_t>(it.row())];
          int gc = global_to_free[static_cast<std::size_t>(it.col())];
          if (gr >= 0 && gc >= 0) {
            triplets.emplace_back(gr, gc, it.value());
          }
        }
      }
      SparseMatrix K_ff(static_cast<int>(num_free), static_cast<int>(num_free));
      K_ff.setFromTriplets(triplets.begin(), triplets.end());
      return K_ff;
    };

    // Scale the convergence tolerance to this step's own force scale,
    // same reasoning as foam::solve_uniaxial_point: the residual is a
    // force (N, or N per unit out-of-plane thickness here), and its
    // natural magnitude depends on the material and mesh size, not a
    // single constant that would suit every case.
    Real reference_force = residual(free_guess).norm();
    NewtonRaphsonSolver::Options options;
    options.tolerance = std::max(reference_force, Real(1.0)) * 1e-10;
    options.max_iterations = 40;
    NewtonRaphsonSolver solver(options);

    auto result = solver.solve(residual, tangent, free_guess);
    if (!result.converged) {
      break;  // report what converged rather than silently extrapolating
    }
    free_guess = result.solution;

    VectorX full_disp = build_full_displacement(t, free_guess);
    VectorX f_int_full = assembler.assemble_residual(full_disp);
    Real total_reaction = 0.0;
    for (auto gdof : loaded_dofs) {
      total_reaction += f_int_full(static_cast<Eigen::Index>(gdof));
    }

    Real strain = std::abs(t * applied_target) / height;
    Real stress = total_reaction / width;  // unit out-of-plane thickness
    curve.emplace_back(strain, stress);
  }

  return curve;
}

DirichletBCs uniaxial_compression_bcs(int nx, int ny, Real top_y_displacement) {
  DirichletBCs bcs;
  for (int i = 0; i <= nx; ++i) {
    std::size_t bottom_id = static_cast<std::size_t>(i);
    std::size_t top_id = static_cast<std::size_t>(ny * (nx + 1) + i);
    bcs.push_back(DirichletBC{bottom_id, 1, 0.0});
    bcs.push_back(DirichletBC{top_id, 1, top_y_displacement});
  }
  bcs.push_back(DirichletBC{0, 0, 0.0});  // pin bottom-left corner's x
  return bcs;
}

}  // namespace foam
