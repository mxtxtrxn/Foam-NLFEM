#pragma once
// Shared linear-algebra aliases. Everything in the project should go
// through these typedefs instead of spelling out Eigen types directly,
// so the underlying scalar/precision or backend can be changed later
// from a single place.

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace foam {

using Real = double;

using Vector2 = Eigen::Matrix<Real, 2, 1>;
using Vector3 = Eigen::Matrix<Real, 3, 1>;
using VectorX = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

using Matrix2 = Eigen::Matrix<Real, 2, 2>;
using Matrix3 = Eigen::Matrix<Real, 3, 3>;
using MatrixX = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

// Voigt-notation tangent moduli (3D): 6x6 symmetric matrix relating
// stress and strain increments in vector form.
using Matrix6 = Eigen::Matrix<Real, 6, 6>;

using SparseMatrix = Eigen::SparseMatrix<Real>;
using Triplet = Eigen::Triplet<Real>;

// Convenience: 3x3 identity, used constantly for F = I + grad(u), etc.
inline Matrix3 identity3() { return Matrix3::Identity(); }

}  // namespace foam
