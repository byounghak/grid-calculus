/// \file
/// \brief 2nd-order accurate periodic divergence on a `Field<Vec3d>`.
///
/// Computes $\nabla \cdot \mathbf{V}$ pointwise using
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

/// \brief Returns $\nabla \cdot \mathbf{V}$ as a fresh, same-shape
///        `Field<double>`.
///
/// 2nd-order central differences:
/// $\partial V_x / \partial x + \partial V_y / \partial y + \partial V_z / \partial z$,
/// each component differentiated along its own axis. Wrapping at the domain
/// boundary is the responsibility of the input `Field`'s `Policy` (default
/// `IndexPolicy::Periodic`); `divergence` itself makes no boundary decisions.
/// \param vfield  Input vector field, components stored as `Vec3d` per grid point.
/// \returns A new `Field<double>` holding $\nabla \cdot \mathbf{V}$ at every grid point.
/// \since 0.2.0
inline core::Field<double> divergence(const core::Field<core::Vec3d>& vfield) {
  using Coeffs = stencil::FirstDerivative<2>;
  const auto& grid = vfield.getGrid();
  const auto& h = grid.getCellSize();
  const double inv_hx = 1.0 / h(0);
  const double inv_hy = 1.0 / h(1);
  const double inv_hz = 1.0 / h(2);

  core::Field<double> result(grid);
  for (int k = 0; k < grid.getNz(); ++k) {
    for (int j = 0; j < grid.getNy(); ++j) {
      for (int i = 0; i < grid.getNx(); ++i) {
        double dvx_dx = 0.0;
        double dvy_dy = 0.0;
        double dvz_dz = 0.0;
        for (std::size_t m = 0; m < Coeffs::offsets.size(); ++m) {
          const int off = Coeffs::offsets[m];
          const double w = Coeffs::weights[m];
          dvx_dx += w * vfield(i + off, j, k)(0);
          dvy_dy += w * vfield(i, j + off, k)(1);
          dvz_dz += w * vfield(i, j, k + off)(2);
        }
        result(i, j, k) = dvx_dx * inv_hx + dvy_dy * inv_hy + dvz_dz * inv_hz;
      }
    }
  }
  return result;
}

}  // namespace gridcalc::diff
