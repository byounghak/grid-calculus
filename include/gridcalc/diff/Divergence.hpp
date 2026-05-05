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
#include <gridcalc/diff/detail/PreconditionAxisExtent.hpp>
#include <gridcalc/stencil/FirstDerivative.hpp>

namespace gridcalc::diff {

/// \brief Returns $\nabla \cdot \mathbf{V}$ as a fresh, same-shape
///        `Field<double>`, using a central-difference stencil of
///        accuracy `Order`.
///
/// Computes $\partial V_x / \partial x + \partial V_y / \partial y + \partial V_z / \partial z$,
/// with each component differentiated along its own axis using the same
/// `stencil::FirstDerivative<Order>` weight table, scaled by the per-axis
/// $1/h_a$. Wrapping at the domain boundary is the responsibility of the
/// input `Field`'s `Policy` (default `IndexPolicy::Periodic`); `divergence`
/// itself makes no boundary decisions.
///
/// `Order` defaults to `2` so callers from Phase 2 onward
/// (`divergence(vfield)`) keep working unchanged. Set explicitly for
/// higher accuracy: `divergence<4>(vfield)` (Phase 7).
/// \tparam Order  Accuracy order of the central-difference stencil. `2`
///                (Phase 2) and `4` (Phase 7) are specialized.
/// \param vfield  Input vector field, components stored as `Vec3d` per grid point.
/// \returns A new `Field<double>` holding $\nabla \cdot \mathbf{V}$ at every grid point.
/// \throws std::invalid_argument if any axis of `vfield`'s grid has
///         extent `N < 2 * FirstDerivative<Order>::radius + 1`. See
///         `diff::detail::requireAxisExtent`.
/// \since 0.2.0 (function); 0.7.0 (`Order` parameter); 0.14.1 (precondition).
template <int Order = 2>
inline core::Field<double> divergence(const core::Field<core::Vec3d>& vfield) {
  using Coeffs = stencil::FirstDerivative<Order>;
  const auto& grid = vfield.getGrid();
  detail::requireAxisExtent("x", grid.getNx(), Coeffs::radius);
  detail::requireAxisExtent("y", grid.getNy(), Coeffs::radius);
  detail::requireAxisExtent("z", grid.getNz(), Coeffs::radius);
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

/// \brief Returns the rank-lowering divergence of a `Field<Mat3d>` as a
///        fresh `Field<Vec3d>`, using a central-difference stencil of
///        accuracy `Order`.
///
/// Index convention: `(div M)_i = ∂_j M(i, j)` (sum over the column / axis
/// index, matching the rank-2 gradient convention `M(i, j) = ∂_j v_i`
/// established by Phase 13's `gradient(Field<Vec3d>)` overload). Under
/// these matched conventions, `divergence(gradient(v))_i = ∂_j ∂_j v_i` is
/// the i-th Cartesian component of the vector Laplacian.
///
/// Implementation: for each output vector component `i`, the stencil is
/// applied along axis `j` to the `M(i, j)` matrix entry; the three axis
/// terms are then summed. Wrapping at the domain boundary is the
/// responsibility of the input `Field`'s `Policy` (default
/// `IndexPolicy::Periodic`); the divergence itself makes no boundary
/// decisions.
/// \tparam Order  Accuracy order of the central-difference stencil. `2`
///                (Phase 2) and `4` (Phase 7) are specialized.
/// \param mfield  Input rank-2 tensor field.
/// \returns A new `Field<Vec3d>` holding `(div M)_i = ∂_j M(i, j)` at every grid point.
/// \throws std::invalid_argument if any axis of `mfield`'s grid has
///         extent `N < 2 * FirstDerivative<Order>::radius + 1`. See
///         `diff::detail::requireAxisExtent`.
/// \since 0.13.0 (function); 0.14.1 (precondition).
template <int Order = 2>
inline core::Field<core::Vec3d> divergence(const core::Field<core::Mat3d>& mfield) {
  using Coeffs = stencil::FirstDerivative<Order>;
  const auto& grid = mfield.getGrid();
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
        core::Vec3d v;
        for (int comp = 0; comp < 3; ++comp) {
          double dx = 0.0;
          double dy = 0.0;
          double dz = 0.0;
          for (std::size_t m = 0; m < Coeffs::offsets.size(); ++m) {
            const int off = Coeffs::offsets[m];
            const double w = Coeffs::weights[m];
            dx += w * mfield(i + off, j, k)(comp, 0);
            dy += w * mfield(i, j + off, k)(comp, 1);
            dz += w * mfield(i, j, k + off)(comp, 2);
          }
          v(comp) = dx * inv_hx + dy * inv_hy + dz * inv_hz;
        }
        result(i, j, k) = v;
      }
    }
  }
  return result;
}

}  // namespace gridcalc::diff
