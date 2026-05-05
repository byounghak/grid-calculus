/// \file
/// \brief Tensor-product helper for multi-axis partial derivatives.
///
/// Implementation detail of the `gridcalc::diff` namespace. Computes
/// $\partial^{N_x+N_y+N_z}\psi / \partial x^{N_x}\partial y^{N_y}\partial z^{N_z}$
/// in a single fused pass over the grid by summing over the outer product of
/// three 1D central-difference weight tables. The Phase 8 named functions
/// (`dxx`, `dxy`, `d3dx2dy`, `d4dx2dy2`, ...) are 2-line wrappers around
/// `mixedDerivative<Nx, Ny, Nz, Order>`.
/// \since 0.8.0

#pragma once

#include <array>
#include <cstddef>

#include <gridcalc/core/Field.hpp>
#include <gridcalc/diff/detail/PreconditionAxisExtent.hpp>
#include <gridcalc/stencil/CentralDifference.hpp>
#include <gridcalc/stencil/FirstDerivative.hpp>
#include <gridcalc/stencil/FourthDerivative.hpp>
#include <gridcalc/stencil/ThirdDerivative.hpp>

namespace gridcalc::diff::detail {

/// \brief Maps a per-axis derivative count `N` and accuracy `Order` to the
///        corresponding 1D central-difference weight table.
///
/// `N = 0` is the identity (single sample, weight 1) so that an axis not
/// involved in the partial contributes a trivial 1-element factor to the
/// outer product. `N = 1..4` aliases the corresponding `gridcalc::stencil`
/// table. The primary template is undefined; instantiating with `N` outside
/// `0..4` is a compile error pointing at the call site.
/// \tparam N      Per-axis derivative count (0..4).
/// \tparam Order  Accuracy order of the underlying central-difference table.
/// \since 0.8.0
template <int N, int Order>
struct AxisStencil;

/// \brief `N = 0` specialization: identity (no derivative on this axis).
/// \since 0.8.0
template <int Order>
struct AxisStencil<0, Order> {
  static constexpr int radius = 0;
  static constexpr std::array<int, 1> offsets = {0};
  static constexpr std::array<double, 1> weights = {1.0};
};

/// \brief `N = 1` specialization: first-derivative table.
/// \since 0.8.0
template <int Order>
struct AxisStencil<1, Order> : stencil::FirstDerivative<Order> {};

/// \brief `N = 2` specialization: second-derivative table.
/// \since 0.8.0
template <int Order>
struct AxisStencil<2, Order> : stencil::Coefficients<Order> {};

/// \brief `N = 3` specialization: third-derivative table.
/// \since 0.8.0
template <int Order>
struct AxisStencil<3, Order> : stencil::ThirdDerivative<Order> {};

/// \brief `N = 4` specialization: fourth-derivative table.
/// \since 0.8.0
template <int Order>
struct AxisStencil<4, Order> : stencil::FourthDerivative<Order> {};

/// \brief Computes the multi-axis partial derivative
///        $\partial^{N_x+N_y+N_z}\psi / \partial x^{N_x}\partial y^{N_y}\partial z^{N_z}$
///        as a single fused pass over the grid.
///
/// At each grid point the helper sums over the outer product of three 1D
/// central-difference weight tables (one per axis). The per-grid-point work
/// is $\prod_a (2 r_a + 1)$ reads where $r_a$ is the radius of the axis-`a`
/// table; an axis with $N_a = 0$ contributes a single read at offset 0.
/// Periodic wrap is delegated to `Field::Policy`.
///
/// Per-axis scaling is $1/h_a^{N_a}$ at the consumer side (the weight
/// tables themselves are unscaled), so anisotropic per-axis spacing is
/// handled without any change to the weight tables.
/// \tparam Nx     X-axis derivative count (0..4).
/// \tparam Ny     Y-axis derivative count (0..4).
/// \tparam Nz     Z-axis derivative count (0..4).
/// \tparam Order  Accuracy order of every per-axis table (must match the
///                supported set on `AxisStencil<N, Order>` for each `N > 0`).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding the requested partial at every
///          grid point.
/// \throws std::invalid_argument if any axis of `psi`'s grid with a
///         non-zero per-axis derivative count has extent
///         `N < 2 * AxisStencil<N_a, Order>::radius + 1`. Axes with
///         per-axis derivative count `0` accept any positive `N` (their
///         stencil radius is `0`). See `diff::detail::requireAxisExtent`.
/// \since 0.8.0 (function); 0.14.1 (precondition).
template <int Nx, int Ny, int Nz, int Order>
inline core::Field<double> mixedDerivative(const core::Field<double>& psi) {
  using Sx = AxisStencil<Nx, Order>;
  using Sy = AxisStencil<Ny, Order>;
  using Sz = AxisStencil<Nz, Order>;

  const auto& grid = psi.getGrid();
  requireAxisExtent("x", grid.getNx(), Sx::radius);
  requireAxisExtent("y", grid.getNy(), Sy::radius);
  requireAxisExtent("z", grid.getNz(), Sz::radius);
  const auto& h = grid.getCellSize();

  double scale = 1.0;
  for (int n = 0; n < Nx; ++n) scale /= h(0);
  for (int n = 0; n < Ny; ++n) scale /= h(1);
  for (int n = 0; n < Nz; ++n) scale /= h(2);

  core::Field<double> result(grid);
  for (int k = 0; k < grid.getNz(); ++k) {
    for (int j = 0; j < grid.getNy(); ++j) {
      for (int i = 0; i < grid.getNx(); ++i) {
        double sum = 0.0;
        for (std::size_t mz = 0; mz < Sz::offsets.size(); ++mz) {
          const int oz = Sz::offsets[mz];
          const double wz = Sz::weights[mz];
          for (std::size_t my = 0; my < Sy::offsets.size(); ++my) {
            const int oy = Sy::offsets[my];
            const double wy = Sy::weights[my];
            const double wyz = wy * wz;
            for (std::size_t mx = 0; mx < Sx::offsets.size(); ++mx) {
              sum += Sx::weights[mx] * wyz *
                     psi(i + Sx::offsets[mx], j + oy, k + oz);
            }
          }
        }
        result(i, j, k) = sum * scale;
      }
    }
  }
  return result;
}

}  // namespace gridcalc::diff::detail
