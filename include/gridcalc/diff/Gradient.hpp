/// \file
/// \brief 2nd-order accurate periodic gradient on a `Field<double>`.
///
/// Computes $\nabla f$ pointwise using
/// `gridcalc::stencil::FirstDerivative<2>` along each axis, scaled by the
/// per-axis $1/h$. Boundary handling is delegated to the input `Field`'s
/// index policy (`Periodic` at Phase 2).
/// \since 0.2.0

#pragma once

#include <cstddef>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/stencil/FirstDerivative.hpp>

namespace gridcalc::diff {

/// \brief Returns $\nabla f$ as a fresh, same-shape `Field<Vec3d>`.
///
/// 2nd-order central differences, applied independently along x, y, and z.
/// At each grid point the result holds the Cartesian triple
/// $(\partial f / \partial x,\ \partial f / \partial y,\ \partial f / \partial z)$.
/// Wrapping at the domain boundary is the responsibility of the input
/// `Field`'s `Policy` (default `IndexPolicy::Periodic`); `gradient` itself
/// makes no boundary decisions.
/// \param field  Input scalar field.
/// \returns A new `Field<Vec3d>` holding $\nabla f$ at every grid point.
/// \since 0.2.0
inline core::Field<core::Vec3d> gradient(const core::Field<double>& field) {
  using Coeffs = stencil::FirstDerivative<2>;
  const auto& grid = field.getGrid();
  const auto& h = grid.getCellSize();
  const double inv_hx = 1.0 / h(0);
  const double inv_hy = 1.0 / h(1);
  const double inv_hz = 1.0 / h(2);

  core::Field<core::Vec3d> result(grid);
  for (int k = 0; k < grid.getNz(); ++k) {
    for (int j = 0; j < grid.getNy(); ++j) {
      for (int i = 0; i < grid.getNx(); ++i) {
        double dx = 0.0;
        double dy = 0.0;
        double dz = 0.0;
        for (std::size_t m = 0; m < Coeffs::offsets.size(); ++m) {
          const int off = Coeffs::offsets[m];
          const double w = Coeffs::weights[m];
          dx += w * field(i + off, j, k);
          dy += w * field(i, j + off, k);
          dz += w * field(i, j, k + off);
        }
        result(i, j, k) = core::Vec3d(dx * inv_hx, dy * inv_hy, dz * inv_hz);
      }
    }
  }
  return result;
}

}  // namespace gridcalc::diff
