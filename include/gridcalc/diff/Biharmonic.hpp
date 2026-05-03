/// \file
/// \brief Biharmonic operator $\nabla^4$ on a periodic `Field<double>`.
///
/// Direct fused-pass implementation (not `laplacian(laplacian(psi))`):
/// $\nabla^4 \psi = \partial_x^4 \psi + \partial_y^4 \psi + \partial_z^4 \psi
///                  + 2(\partial_x^2 \partial_y^2 \psi
///                      + \partial_x^2 \partial_z^2 \psi
///                      + \partial_y^2 \partial_z^2 \psi)$.
/// The composition `laplacian(laplacian(psi))` would inherit the stride-2
/// leading-error behavior already documented for `divergence(gradient(psi))`
/// — same convergence order but a different leading constant. The direct
/// stencil sidesteps that drift and lands the roadmap acceptance bar
/// (biharmonic of $\sin(\mathbf{k}\cdot\mathbf{x})$ recovers $|\mathbf{k}|^4$)
/// at the published constant.
/// \since 0.8.0

#pragma once

#include <cstddef>

#include <gridcalc/core/Field.hpp>
#include <gridcalc/stencil/CentralDifference.hpp>
#include <gridcalc/stencil/FourthDerivative.hpp>

namespace gridcalc::diff {

/// \brief Returns $\nabla^4 \psi$ as a fresh same-shape field, using a
///        direct six-term tensor-product stencil.
///
/// Sums the three single-axis 4th derivatives ($\partial_x^4$, $\partial_y^4$,
/// $\partial_z^4$) and twice the three mixed (2, 2) cross terms
/// ($\partial_x^2 \partial_y^2$, $\partial_x^2 \partial_z^2$,
/// $\partial_y^2 \partial_z^2$) in a single pass over the grid. Per-axis
/// scaling is $1/h_a^4$ for the single-axis terms and $1/(h_a^2 h_b^2)$ for
/// the mixed terms — anisotropic per-axis spacing handled without weight-table
/// modification. Periodic wrap is delegated to `Field::Policy`.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\nabla^4 \psi$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> biharmonic(const core::Field<double>& psi) {
  using S2 = stencil::Coefficients<Order>;
  using S4 = stencil::FourthDerivative<Order>;

  const auto& grid = psi.getGrid();
  const auto& h = grid.getCellSize();
  const double inv_hx2 = 1.0 / (h(0) * h(0));
  const double inv_hy2 = 1.0 / (h(1) * h(1));
  const double inv_hz2 = 1.0 / (h(2) * h(2));
  const double inv_hx4 = inv_hx2 * inv_hx2;
  const double inv_hy4 = inv_hy2 * inv_hy2;
  const double inv_hz4 = inv_hz2 * inv_hz2;

  core::Field<double> result(grid);
  for (int k = 0; k < grid.getNz(); ++k) {
    for (int j = 0; j < grid.getNy(); ++j) {
      for (int i = 0; i < grid.getNx(); ++i) {
        double dx4 = 0.0;
        double dy4 = 0.0;
        double dz4 = 0.0;
        for (std::size_t m = 0; m < S4::offsets.size(); ++m) {
          const int off = S4::offsets[m];
          const double w = S4::weights[m];
          dx4 += w * psi(i + off, j, k);
          dy4 += w * psi(i, j + off, k);
          dz4 += w * psi(i, j, k + off);
        }

        double dx2y2 = 0.0;
        double dx2z2 = 0.0;
        double dy2z2 = 0.0;
        for (std::size_t mb = 0; mb < S2::offsets.size(); ++mb) {
          const int ob = S2::offsets[mb];
          const double wb = S2::weights[mb];
          for (std::size_t ma = 0; ma < S2::offsets.size(); ++ma) {
            const int oa = S2::offsets[ma];
            const double w_ab = S2::weights[ma] * wb;
            dx2y2 += w_ab * psi(i + oa, j + ob, k);
            dx2z2 += w_ab * psi(i + oa, j, k + ob);
            dy2z2 += w_ab * psi(i, j + oa, k + ob);
          }
        }

        result(i, j, k) =
            dx4 * inv_hx4 + dy4 * inv_hy4 + dz4 * inv_hz4 +
            2.0 * (dx2y2 * inv_hx2 * inv_hy2 + dx2z2 * inv_hx2 * inv_hz2 +
                   dy2z2 * inv_hy2 * inv_hz2);
      }
    }
  }
  return result;
}

/// \brief Alias `d4<Order>(psi)` for `biharmonic<Order>(psi)`.
///
/// Completes the d-prefix family at the contracted-rank level: with no axis
/// suffix, `d<rank>` denotes the contracted scalar operator on an even rank
/// (`d4` = $\nabla^4$). `laplacian` is *not* renamed and no `d2` alias is
/// added at Phase 8 — Phase 1's API stays put.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\nabla^4 \psi$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d4(const core::Field<double>& psi) {
  return biharmonic<Order>(psi);
}

}  // namespace gridcalc::diff
