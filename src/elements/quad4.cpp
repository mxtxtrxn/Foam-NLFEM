#include "foam/elements/quad4.hpp"

#include <cmath>
#include <stdexcept>

#include "foam/kinematics/deformation_gradient.hpp"
#include "foam/materials/material_base.hpp"

namespace foam {

namespace {

// Standard 2x2 Gauss quadrature points for the [-1,1]x[-1,1] reference
// square; each carries weight 1 (sum of weights = 4 = reference area).
std::array<std::array<Real, 2>, 4> gauss_points_2x2() {
  Real gp = 1.0 / std::sqrt(3.0);
  return {{{-gp, -gp}, {gp, -gp}, {gp, gp}, {-gp, gp}}};
}

// Bilinear Quad4 shape-function derivatives w.r.t. natural coordinates
// (xi, eta), node order bottom-left/bottom-right/top-right/top-left.
Eigen::Matrix<Real, 4, 2> shape_gradients_natural(Real xi, Real eta) {
  Eigen::Matrix<Real, 4, 2> dN;
  dN(0, 0) = -(1 - eta) / 4.0;
  dN(0, 1) = -(1 - xi) / 4.0;
  dN(1, 0) = (1 - eta) / 4.0;
  dN(1, 1) = -(1 + xi) / 4.0;
  dN(2, 0) = (1 + eta) / 4.0;
  dN(2, 1) = (1 + xi) / 4.0;
  dN(3, 0) = -(1 + eta) / 4.0;
  dN(3, 1) = (1 - xi) / 4.0;
  return dN;
}

}  // namespace

std::array<Quad4Element::GaussPointData, 4> Quad4Element::compute_gauss_point_data(
    const MatrixX& reference_positions, const MatrixX& nodal_displacements) const {
  if (reference_positions.rows() != 4 || reference_positions.cols() != 2 ||
      nodal_displacements.rows() != 4 || nodal_displacements.cols() != 2) {
    throw std::invalid_argument(
        "Quad4Element: reference_positions and nodal_displacements must each be 4x2 (4 nodes, "
        "2 dofs/node)");
  }

  std::array<GaussPointData, 4> gauss_data;
  auto gauss_points = gauss_points_2x2();

  for (int g = 0; g < 4; ++g) {
    Real xi = gauss_points[g][0];
    Real eta = gauss_points[g][1];
    Eigen::Matrix<Real, 4, 2> dN_dxi = shape_gradients_natural(xi, eta);

    // Reference Jacobian J0(i,j) = sum_a X_a(i) * dN_a/dxi_j, mapping
    // dX = J0 * d(xi,eta).
    Matrix2 J0 = Matrix2::Zero();
    for (int a = 0; a < 4; ++a) {
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          J0(i, j) += reference_positions(a, i) * dN_dxi(a, j);
        }
      }
    }

    Real detJ0 = J0.determinant();
    if (detJ0 <= 0.0) {
      throw std::runtime_error(
          "Quad4Element: non-positive reference Jacobian determinant (degenerate element, or "
          "nodes not ordered bottom-left/bottom-right/top-right/top-left CCW)");
    }
    Matrix2 J0_inv_T = J0.inverse().transpose();

    // Reference shape-function gradients dN/dX = J0^-T * dN/dxi.
    GaussPointData data;
    data.detJ0 = detJ0;
    for (int a = 0; a < 4; ++a) {
      Eigen::Matrix<Real, 2, 1> grad_xi = dN_dxi.row(a).transpose();
      data.dN_dX.row(a) = (J0_inv_T * grad_xi).transpose();
    }

    data.F = compute_deformation_gradient(data.dN_dX, nodal_displacements);

    Real F11 = data.F(0, 0), F12 = data.F(0, 1), F21 = data.F(1, 0), F22 = data.F(1, 1);
    data.B_L.setZero();
    for (int a = 0; a < 4; ++a) {
      Real dNa1 = data.dN_dX(a, 0);
      Real dNa2 = data.dN_dX(a, 1);
      // x-dof column (node a's x-displacement, index 2a)
      data.B_L(0, 2 * a) = F11 * dNa1;
      data.B_L(1, 2 * a) = F12 * dNa2;
      data.B_L(2, 2 * a) = F11 * dNa2 + F12 * dNa1;
      // y-dof column (node a's y-displacement, index 2a+1)
      data.B_L(0, 2 * a + 1) = F21 * dNa1;
      data.B_L(1, 2 * a + 1) = F22 * dNa2;
      data.B_L(2, 2 * a + 1) = F21 * dNa2 + F22 * dNa1;
    }

    gauss_data[g] = data;
  }

  return gauss_data;
}

VectorX Quad4Element::internal_force(const MatrixX& reference_positions,
                                      const MatrixX& nodal_displacements,
                                      const MaterialModel& material) const {
  auto gauss_data = compute_gauss_point_data(reference_positions, nodal_displacements);

  VectorX f = VectorX::Zero(8);
  for (const auto& gp : gauss_data) {
    Matrix3 S_full = second_piola_kirchhoff_stress(material, gp.F);
    Eigen::Matrix<Real, 3, 1> S_reduced;
    S_reduced << S_full(0, 0), S_full(1, 1), S_full(0, 1);
    // Gauss weight is 1 for this 2x2 rule, so just multiply by detJ0.
    f += gp.B_L.transpose() * S_reduced * gp.detJ0;
  }
  return f;
}

MatrixX Quad4Element::tangent_stiffness(const MatrixX& reference_positions,
                                         const MatrixX& nodal_displacements,
                                         const MaterialModel& material) const {
  auto gauss_data = compute_gauss_point_data(reference_positions, nodal_displacements);

  MatrixX K = MatrixX::Zero(8, 8);
  for (const auto& gp : gauss_data) {
    Matrix3 S_full = second_piola_kirchhoff_stress(material, gp.F);
    Matrix6 D_full = material.tangent_moduli(gp.F);

    // Material part: reduce D (6x6, full 3D Voigt) to the 3x3 in-plane
    // block at Voigt indices {11,22,12} = {0,1,3}. This is valid
    // because E33 = E13 = E23 = 0 IDENTICALLY for plane strain (no
    // out-of-plane displacement, no z-dependence), not just small: with
    // dE33 = dE13 = dE23 = 0 always, those columns/rows of D never
    // contribute, regardless of their actual values.
    static constexpr int kIdx[3] = {0, 1, 3};
    Eigen::Matrix<Real, 3, 3> D_reduced;
    for (int p = 0; p < 3; ++p) {
      for (int q = 0; q < 3; ++q) {
        D_reduced(p, q) = D_full(kIdx[p], kIdx[q]);
      }
    }

    K += gp.B_L.transpose() * D_reduced * gp.B_L * gp.detJ0;

    // Geometric (initial-stress) stiffness: a SCALAR g_ab = grad_a^T *
    // S_inplane * grad_b added to the (a,x)-(b,x) and (a,y)-(b,y)
    // entries only (zero cross-coupling between x and y), per the
    // standard Total Lagrangian derivation (see docs/roadmap.md note
    // for this step, or derive it directly: d(F*S*dN_a)/d(u_b,l) =
    // delta_{kl} * (dN_a . S . dN_b) since dF_kj/d(u_b,l) = delta_kl *
    // dN_b[j]).
    Matrix2 S_inplane;
    S_inplane << S_full(0, 0), S_full(0, 1), S_full(0, 1), S_full(1, 1);
    for (int a = 0; a < 4; ++a) {
      Eigen::Matrix<Real, 2, 1> grad_a = gp.dN_dX.row(a).transpose();
      for (int b = 0; b < 4; ++b) {
        Eigen::Matrix<Real, 2, 1> grad_b = gp.dN_dX.row(b).transpose();
        Real g_ab = (grad_a.transpose() * S_inplane * grad_b)(0, 0) * gp.detJ0;
        K(2 * a, 2 * b) += g_ab;
        K(2 * a + 1, 2 * b + 1) += g_ab;
      }
    }
  }
  return K;
}

}  // namespace foam
