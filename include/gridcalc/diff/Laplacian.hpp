/// \file
/// \brief 2nd-order accurate periodic Laplacian on a `Field<double>`.
///
/// Computes $\nabla^2 f$ pointwise using
/// `gridcalc::stencil::Coefficients<2>` along each axis, scaled by the
/// per-axis $1/h^2$. Boundary handling is delegated to the input `Field`'s
/// index policy (`Periodic` at Phase 1).
/// \since 0.1.0

#pragma once

#include <cstddef>

#include <gridcalc/core/Field.hpp>
#include <gridcalc/stencil/CentralDifference.hpp>

namespace gridcalc::diff {

/// \brief Returns $\nabla^2 f$ as a fresh, same-shape `Field<double>`,
///        using a central-difference stencil of accuracy `Order`.
///
/// Each Cartesian axis is treated independently with the same
/// `stencil::Coefficients<Order>` weight table, scaled by that axis's
/// `1/h_a^2`; the three per-axis sums are added at each grid point.
/// Wrapping at the domain boundary is the responsibility of the input
/// `Field`'s `Policy` (default `IndexPolicy::Periodic`); `laplacian`
/// itself makes no boundary decisions.
///
/// `Order` defaults to `2` so callers from Phase 1 onward (`laplacian(field)`)
/// keep working unchanged. Set explicitly for higher accuracy:
/// `laplacian<4>(field)` (Phase 7). Unsupported orders trip a compile
/// error inside `Coefficients<Order>`.
/// \tparam Order  Accuracy order of the central-difference stencil. Must
///                be a supported even integer; `2` (Phase 1) and `4`
///                (Phase 7) are specialized.
/// \param field   Input scalar field.
/// \returns A new `Field<double>` holding $\nabla^2 f$ at every grid point.
/// \since 0.1.0 (function); 0.7.0 (`Order` parameter).
template <int Order = 2>
inline core::Field<double> laplacian(const core::Field<double>& field) {
  using Coeffs = stencil::Coefficients<Order>;
  const auto& grid = field.getGrid();
  const auto& h = grid.getCellSize();
  const double inv_hx2 = 1.0 / (h(0) * h(0));
  const double inv_hy2 = 1.0 / (h(1) * h(1));
  const double inv_hz2 = 1.0 / (h(2) * h(2));

  core::Field<double> result(grid);
  for (int k = 0; k < grid.getNz(); ++k) {
    for (int j = 0; j < grid.getNy(); ++j) {
      for (int i = 0; i < grid.getNx(); ++i) {
        double dxx = 0.0;
        double dyy = 0.0;
        double dzz = 0.0;
        for (std::size_t m = 0; m < Coeffs::offsets.size(); ++m) {
          const int off = Coeffs::offsets[m];
          const double w = Coeffs::weights[m];
          dxx += w * field(i + off, j, k);
          dyy += w * field(i, j + off, k);
          dzz += w * field(i, j, k + off);
        }
        result(i, j, k) = dxx * inv_hx2 + dyy * inv_hy2 + dzz * inv_hz2;
      }
    }
  }
  return result;
}

}  // namespace gridcalc::diff
