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
#include <gridcalc/diff/detail/PreconditionAxisExtent.hpp>
#include <gridcalc/stencil/FirstDerivative.hpp>

namespace gridcalc::diff {

/// \brief Returns $\nabla f$ as a fresh, same-shape `Field<Vec3d>`,
///        using a central-difference stencil of accuracy `Order`.
///
/// At each grid point the result holds the Cartesian triple
/// $(\partial f / \partial x,\ \partial f / \partial y,\ \partial f / \partial z)$.
/// Each component uses the same `stencil::FirstDerivative<Order>` weight
/// table along its own axis, scaled by the per-axis $1/h_a$. Wrapping at
/// the domain boundary is the responsibility of the input `Field`'s
/// `Policy` (default `IndexPolicy::Periodic`); `gradient` itself makes
/// no boundary decisions.
///
/// `Order` defaults to `2` so callers from Phase 2 onward
/// (`gradient(field)`) keep working unchanged. Set explicitly for higher
/// accuracy: `gradient<4>(field)` (Phase 7).
/// \tparam Order  Accuracy order of the central-difference stencil. `2`
///                (Phase 2) and `4` (Phase 7) are specialized.
/// \param field  Input scalar field.
/// \returns A new `Field<Vec3d>` holding $\nabla f$ at every grid point.
/// \throws std::invalid_argument if any axis of `field`'s grid has
///         extent `N < 2 * FirstDerivative<Order>::radius + 1`. See
///         `diff::detail::requireAxisExtent`.
/// \since 0.2.0 (function); 0.7.0 (`Order` parameter); 0.14.1 (precondition).
template <int Order = 2>
inline core::Field<core::Vec3d> gradient(const core::Field<double>& field) {
  using Coeffs = stencil::FirstDerivative<Order>;
  const auto& grid = field.getGrid();
  detail::requireAxisExtent("x", grid.getNx(), Coeffs::radius);
  detail::requireAxisExtent("y", grid.getNy(), Coeffs::radius);
  detail::requireAxisExtent("z", grid.getNz(), Coeffs::radius);
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

/// \brief Returns the rank-2 gradient of a `Field<Vec3d>` as a fresh
///        `Field<Mat3d>`, using a central-difference stencil of accuracy `Order`.
///
/// Index convention: `M(i, j) = ∂_j v_i`. Row `i` selects the vector
/// component; column `j` selects the derivative axis. Under this
/// convention `trace(M) = ∂_x v_x + ∂_y v_y + ∂_z v_z = ∇·v`, and
/// `divergence(M)_i = ∂_j M(i, j)` gives the standard rank-lowering
/// divergence used by Phase 13's `divergence(Field<Mat3d>)` overload.
///
/// Implementation is a direct extension of the scalar gradient: for each
/// vector component `c ∈ {x, y, z}`, the central-difference stencil is
/// applied along each Cartesian axis and the resulting derivatives populate
/// row `c` of the output matrix. Wrapping at the domain boundary is the
/// responsibility of the input `Field`'s `Policy` (default
/// `IndexPolicy::Periodic`); the gradient itself makes no boundary
/// decisions.
/// \tparam Order  Accuracy order of the central-difference stencil. `2`
///                (Phase 2) and `4` (Phase 7) are specialized.
/// \param field   Input vector field.
/// \returns A new `Field<Mat3d>` holding `∂_j v_i` at every grid point.
/// \throws std::invalid_argument if any axis of `field`'s grid has
///         extent `N < 2 * FirstDerivative<Order>::radius + 1`. See
///         `diff::detail::requireAxisExtent`.
/// \since 0.13.0 (function); 0.14.1 (precondition).
template <int Order = 2>
inline core::Field<core::Mat3d> gradient(const core::Field<core::Vec3d>& field) {
  using Coeffs = stencil::FirstDerivative<Order>;
  const auto& grid = field.getGrid();
  detail::requireAxisExtent("x", grid.getNx(), Coeffs::radius);
  detail::requireAxisExtent("y", grid.getNy(), Coeffs::radius);
  detail::requireAxisExtent("z", grid.getNz(), Coeffs::radius);
  const auto& h = grid.getCellSize();
  const double inv_hx = 1.0 / h(0);
  const double inv_hy = 1.0 / h(1);
  const double inv_hz = 1.0 / h(2);

  core::Field<core::Mat3d> result(grid);
  for (int k = 0; k < grid.getNz(); ++k) {
    for (int j = 0; j < grid.getNy(); ++j) {
      for (int i = 0; i < grid.getNx(); ++i) {
        core::Mat3d M;
        for (int comp = 0; comp < 3; ++comp) {
          double dx = 0.0;
          double dy = 0.0;
          double dz = 0.0;
          for (std::size_t m = 0; m < Coeffs::offsets.size(); ++m) {
            const int off = Coeffs::offsets[m];
            const double w = Coeffs::weights[m];
            dx += w * field(i + off, j, k)(comp);
            dy += w * field(i, j + off, k)(comp);
            dz += w * field(i, j, k + off)(comp);
          }
          M(comp, 0) = dx * inv_hx;
          M(comp, 1) = dy * inv_hy;
          M(comp, 2) = dz * inv_hz;
        }
        result(i, j, k) = M;
      }
    }
  }
  return result;
}

}  // namespace gridcalc::diff
