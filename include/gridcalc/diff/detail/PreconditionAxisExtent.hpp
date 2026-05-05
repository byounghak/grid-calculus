/// \file
/// \brief Per-axis stencil-radius precondition shared by every
///        `gridcalc::diff::*` and `gridcalc::fvm::*` operator that reads
///        off-axis samples through the input `Field`'s periodic index
///        policy.
///
/// Centered finite-difference stencils of radius `r` reach offsets
/// `[-r, r]` along each active axis. On a periodic axis with extent
/// `N <= 2 * r`, the wrap aliases two or more offsets to the same grid
/// sample, silently degrading the stencil away from its tabulated
/// truncation order. This helper lets each operator's entry point
/// reject those grids loudly with a precise message, instead of
/// returning silently-wrong derivatives.
/// \since 0.14.1
///
/// **Implementation detail.** Lives under `diff::detail` so the helper
/// stays out of the public API surface; `fvm::cellLaplacian` reuses it
/// via a relative include.

#pragma once

#include <stdexcept>
#include <string>

namespace gridcalc::diff::detail {

/// \brief Throws `std::invalid_argument` if a periodic axis is in the
///        ambiguous "partial-alias" extent range for the operator's
///        stencil radius.
///
/// On a periodic axis with extent `N`, a centered stencil of radius
/// `radius` reads offsets `[-radius, radius]`. Three regimes:
///
///  - **`N >= 2 * radius + 1`** — every offset lands on a distinct
///    sample after the periodic wrap; the tabulated stencil is exact
///    as advertised.
///  - **`N == 1`** — every offset wraps to index `0`. All `2*radius+1`
///    reads return the single cell value `psi(0)`, and the stencil
///    weights sum to zero by construction (the n-th derivative of a
///    constant is zero), so the result is exactly `0`. This is the
///    *correct* value for a degenerate single-cell axis where there is
///    no spatial variation, and is therefore accepted.
///  - **`2 <= N <= 2 * radius`** — *some* but not all offsets alias.
///    The aliased weights add at the duplicated index; the resulting
///    operator is no longer the tabulated stencil and silently drifts
///    away from the advertised truncation order. This range is
///    rejected.
///
/// When the precondition fails, the message names the axis label, its
/// `N`, the stencil radius, and the required minimum so the caller can
/// fix the grid (or pick a lower-`Order` stencil) without reading
/// source.
///
/// `radius == 0` is a no-op: an axis the operator does not differentiate
/// (e.g., the z axis in `d2dxy<Order>`) has stencil radius `0` through
/// `AxisStencil<0, Order>` and accepts any positive `N`.
///
/// \param axis_label  Short axis tag for the error message ("x", "y", "z").
/// \param N           Periodic axis extent (must be `> 0`; the `Grid`
///                    constructor already rejects non-positive extents).
/// \param radius      Stencil half-width along this axis.
/// \throws std::invalid_argument if `2 <= N <= 2 * radius`.
/// \since 0.14.1
inline void requireAxisExtent(const char* axis_label, int N, int radius) {
  if (radius == 0) return;
  if (N == 1) return;
  const int min_N = 2 * radius + 1;
  if (N < min_N) {
    throw std::invalid_argument(
        "gridcalc: axis " + std::string(axis_label) + " has N=" +
        std::to_string(N) + " but the requested stencil has radius " +
        std::to_string(radius) +
        ", which requires N>=" + std::to_string(min_N) +
        " on a periodic grid (or N=1 for a degenerate axis); "
        "intermediate extents alias the [-r, r] offsets under the "
        "periodic wrap and the stencil silently degrades. Use a larger "
        "grid or a lower-Order stencil.");
  }
}

}  // namespace gridcalc::diff::detail
