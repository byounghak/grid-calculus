/// \file
/// \brief Pointwise contraction primitives for `Field<Mat3d>`.
///
/// Phase 13 ships fully-materialized contractions: every operation
/// allocates a fresh `Field<...>` for its output. Phase 14 may replace
/// these with expression-template / fused-loop variants alongside the
/// `func::evaluate`-over-tensor-fields surface that motivates them; for
/// now, with rank-2 inputs only, the rank-6 functional-eval explosion
/// the ET layer is meant to fend off does not yet appear.
/// \since 0.13.0

#pragma once

#include <cstddef>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/tensor/detail/RequireSameGrid.hpp>

namespace gridcalc::tensor {

/// \brief Returns the pointwise trace of a `Field<Mat3d>` as a fresh
///        `Field<double>`.
///
/// At each grid point: `trace(M)(i, j, k) = M(i, j, k).trace() =
/// M(0, 0) + M(1, 1) + M(2, 2)`. Under the Phase 13 rank-2 gradient
/// convention `M(i, j) = ∂_j v_i`, the trace gives the divergence of the
/// underlying vector field --- the second roadmap-acceptance identity.
/// \param mfield  Input rank-2 tensor field.
/// \returns A new `Field<double>` holding the per-point trace.
/// \since 0.13.0
inline core::Field<double> trace(const core::Field<core::Mat3d>& mfield) {
  const auto& grid = mfield.getGrid();
  core::Field<double> result(grid);
  const std::size_t n = grid.getSize();
  const core::Mat3d* in = mfield.data();
  double* out = result.data();
  for (std::size_t idx = 0; idx < n; ++idx) {
    out[idx] = in[idx].trace();
  }
  return result;
}

/// \brief Returns the pointwise single contraction (matrix product) of
///        two `Field<Mat3d>` operands as a fresh `Field<Mat3d>`.
///
/// At each grid point: `(A · B)(i, j) = A(i, k) B(k, j)` (sum over `k`),
/// i.e., `out(i, j, k) = a(i, j, k) * b(i, j, k)` interpreting the `*`
/// as `Eigen::Matrix3d::operator*`. The two inputs must share the same
/// `Grid` (extents and cell size); a mismatch throws
/// `std::invalid_argument` (added 0.15.0; before that the loop was a
/// silent flat-index OOB / mis-pairing on non-matching grids).
/// \param a  Left operand.
/// \param b  Right operand.
/// \returns A new `Field<Mat3d>` holding the per-point matrix product.
/// \throws std::invalid_argument if `a.getGrid() != b.getGrid()`.
/// \since 0.13.0 (function); 0.15.0 (Grid-equality precondition).
inline core::Field<core::Mat3d> singleContract(const core::Field<core::Mat3d>& a,
                                               const core::Field<core::Mat3d>& b) {
  detail::requireSameGrid(a, b, "tensor::singleContract");
  const auto& grid = a.getGrid();
  core::Field<core::Mat3d> result(grid);
  const std::size_t n = grid.getSize();
  const core::Mat3d* pa = a.data();
  const core::Mat3d* pb = b.data();
  core::Mat3d* po = result.data();
  for (std::size_t idx = 0; idx < n; ++idx) {
    po[idx] = pa[idx] * pb[idx];
  }
  return result;
}

/// \brief Returns the pointwise full (double) contraction of two
///        `Field<Mat3d>` operands as a fresh `Field<double>`.
///
/// At each grid point: `A : B = A(i, j) B(i, j)` (sum over both
/// indices), equivalently `(A.array() * B.array()).sum()` or
/// `trace(Aᵀ B)`. The standard pairing for elastic energy density
/// `½ σ : ε` and similar rank-2 / rank-2 contractions. The two inputs
/// must share the same `Grid`; a mismatch throws
/// `std::invalid_argument` (added 0.15.0; before that the loop was a
/// silent flat-index OOB / mis-pairing on non-matching grids).
/// \param a  Left operand.
/// \param b  Right operand.
/// \returns A new `Field<double>` holding the per-point full contraction.
/// \throws std::invalid_argument if `a.getGrid() != b.getGrid()`.
/// \since 0.13.0 (function); 0.15.0 (Grid-equality precondition).
inline core::Field<double> doubleContract(const core::Field<core::Mat3d>& a,
                                          const core::Field<core::Mat3d>& b) {
  detail::requireSameGrid(a, b, "tensor::doubleContract");
  const auto& grid = a.getGrid();
  core::Field<double> result(grid);
  const std::size_t n = grid.getSize();
  const core::Mat3d* pa = a.data();
  const core::Mat3d* pb = b.data();
  double* po = result.data();
  for (std::size_t idx = 0; idx < n; ++idx) {
    po[idx] = (pa[idx].array() * pb[idx].array()).sum();
  }
  return result;
}

}  // namespace gridcalc::tensor
