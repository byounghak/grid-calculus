/// \file
/// \brief Project-local Eigen typedefs shared across `gridcalc` headers.
/// \since 0.2.0

#pragma once

#include <Eigen/Core>

namespace gridcalc::core {

/// \brief Project-local alias for the 3-component Cartesian vector type.
///
/// Originally declared in `core/Grid.hpp` (Phase 1) and relocated here in
/// Phase 2 once `diff/Gradient.hpp` and `diff/Divergence.hpp` also needed it.
/// The fully-qualified name `gridcalc::core::Vec3d` is unchanged from
/// Phase 1.
/// \since 0.1.0
using Vec3d = Eigen::Vector3d;

/// \brief Project-local alias for the rank-2 3x3 Cartesian tensor type.
///
/// Backed directly by `Eigen::Matrix3d` (standard Eigen module — not the
/// unsupported tensor module). Used by Phase 13's rank-raising
/// `diff::gradient(Field<Vec3d>) -> Field<Mat3d>` (with the convention
/// `M(i, j) = ∂_j v_i`, so the trace gives the divergence) and
/// `diff::divergence(Field<Mat3d>) -> Field<Vec3d>`, and by the
/// contraction primitives in `tensor/Contraction.hpp`. Per Eigen's
/// contract, the default constructor leaves coefficients uninitialized;
/// pass `Mat3d::Zero()` (or any other broadcast value) to the `Field`
/// value-broadcast constructor when zero initialization is needed.
/// \since 0.13.0
using Mat3d = Eigen::Matrix3d;

}  // namespace gridcalc::core
