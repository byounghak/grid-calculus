/// \file
/// \brief Periodic Cartesian-grid integration of a scalar `Field<double>`
///        or a `double`-valued contraction-ET expression.
///
/// Implements the periodic midpoint quadrature rule
/// $I_h = h_x h_y h_z \sum_{i,j,k} f_{ijk}$, which on the equispaced periodic
/// grid is identical to the trapezoidal rule and is spectrally accurate for
/// analytic periodic integrands. Two reduction strategies are exposed via
/// tag dispatch: pairwise recursive summation (default) and Kahan/Neumaier
/// compensated summation.
///
/// Phase 15 adds a fused overload accepting any contraction-ET node from
/// `gridcalc/tensor/Expressions.hpp` whose `value_type` is `double`,
/// avoiding `Field<core::Mat3d>` intermediates when reducing tensor
/// expressions like the linear elastic energy
/// $\int \tfrac12 \boldsymbol{\sigma}:\boldsymbol{\varepsilon}\,dV$.
/// \since 0.3.0 (scalar `Field<double>` overloads); 0.15.0 (ET overload).

#pragma once

#include <cmath>
#include <cstddef>
#include <type_traits>

#include <gridcalc/core/Field.hpp>
#include <gridcalc/tensor/Expressions.hpp>

namespace gridcalc::func {

/// \brief Tag selecting pairwise recursive summation as the reduction.
///
/// Default reducer for `integrate()`. Rounding error grows as
/// $O(\varepsilon \log n)$ where $\varepsilon$ is `DBL_EPSILON` and $n$ is
/// the number of grid points. Cheap and easy to keep deterministic across
/// future parallelization (Phase 20).
/// \since 0.3.0
struct Pairwise {};

/// \brief Tag selecting Kahan/Neumaier compensated summation as the reduction.
///
/// Single-pass per-element compensation; rounding error is independent of
/// $n$ to first order. Cost is ~3-4x the pairwise reducer per element. Use
/// when the field contains values whose magnitudes vary by many orders of
/// magnitude.
/// \since 0.3.0
struct Kahan {};

namespace detail {

/// \brief Pairwise recursive sum of `n` doubles starting at `data`.
///
/// Recurses until the block size is below `kBaseCase`, at which point a
/// naive accumulator runs over the block. For `n` smaller than `kBaseCase`
/// this is equivalent to a naive sum.
/// \param data  Pointer to a contiguous block of `n` doubles.
/// \param n     Number of elements to sum.
/// \returns The pairwise-reduced sum.
/// \since 0.3.0
inline double pairwiseSum(const double* data, std::size_t n) noexcept {
  static constexpr std::size_t kBaseCase = 64;
  if (n <= kBaseCase) {
    double s = 0.0;
    for (std::size_t i = 0; i < n; ++i) s += data[i];
    return s;
  }
  const std::size_t mid = n / 2;
  return pairwiseSum(data, mid) + pairwiseSum(data + mid, n - mid);
}

/// \brief Neumaier-style compensated sum of `n` doubles starting at `data`.
///
/// Variant of Kahan summation that handles the `|sum| < |x|` case by
/// compensating against `x` rather than against `sum`. The compensation
/// term is added back into the result at the end.
/// \param data  Pointer to a contiguous block of `n` doubles.
/// \param n     Number of elements to sum.
/// \returns The compensated sum.
/// \since 0.3.0
inline double neumaierSum(const double* data, std::size_t n) noexcept {
  double sum = 0.0;
  double c = 0.0;  // running compensation
  for (std::size_t i = 0; i < n; ++i) {
    const double x = data[i];
    const double t = sum + x;
    if (std::abs(sum) >= std::abs(x)) {
      c += (sum - t) + x;
    } else {
      c += (x - t) + sum;
    }
    sum = t;
  }
  return sum + c;
}

}  // namespace detail

/// \brief Returns $\int f\, dV$ using pairwise recursive summation (default).
///
/// Computes the discrete periodic-midpoint integral
/// $I_h = h_x h_y h_z \sum_{i,j,k} f_{ijk}$ over all grid points. The
/// pairwise reducer has rounding-error bound $O(\varepsilon \log n)$ vs the
/// naive $O(\varepsilon n)$. At Phase 3 the implementation is
/// single-threaded; the result is therefore trivially deterministic across
/// thread counts. Phase 20 parallelizes the reduction while preserving
/// bit-identical output via fixed-size chunking.
/// \param field  Input scalar field.
/// \param tag    Reduction-strategy tag (defaulted to `Pairwise{}`).
/// \returns The discrete integral as a `double`.
/// \since 0.3.0
inline double integrate(const core::Field<double>& field, Pairwise tag = {}) noexcept {
  (void)tag;
  return field.getGrid().getCellVolume() *
         detail::pairwiseSum(field.data(), field.getSize());
}

/// \brief Returns $\int f\, dV$ using Kahan/Neumaier compensated summation.
///
/// Single-pass per-element compensation; ~bit-perfect accuracy at ~3-4x the
/// cost per element of the pairwise reducer. Use when the field contains
/// values whose magnitudes vary by many orders of magnitude.
/// \param field  Input scalar field.
/// \param tag    `Kahan{}` (selects this overload).
/// \returns The discrete integral as a `double`.
/// \since 0.3.0
inline double integrate(const core::Field<double>& field, Kahan tag) noexcept {
  (void)tag;
  return field.getGrid().getCellVolume() *
         detail::neumaierSum(field.data(), field.getSize());
}

/// \brief Returns $\int \mathrm{expr}\, dV$ for a `double`-valued contraction-ET expression.
///
/// Materializes the per-cell scalar integrand once into a temporary
/// `core::Field<double>` and forwards to the matching scalar overload.
/// The point of the overload is to **avoid** materializing rank-2
/// intermediates (`Field<core::Mat3d>` for $\boldsymbol{\sigma}$,
/// $\boldsymbol{\varepsilon}$, etc.) along the way: the expression's
/// `evalAt(i, j, k)` evaluator pulls the pointwise scalar directly from
/// the Mat3d arithmetic at each cell.
///
/// Equivalent in numerical output to
/// `func::integrate(tensor::expr::materialize(expr), tag)`. The
/// expression's `value_type` must be `double`, enforced via
/// `std::enable_if`. Supported tags are `Pairwise` (default) and `Kahan`.
/// \tparam E    A contraction-ET node type (`tensor::expr::is_expr_v<E>` must be `true`)
///              with `value_type == double`.
/// \tparam Tag  Reduction-strategy tag --- `Pairwise` or `Kahan`.
/// \param expr  The scalar-valued expression to reduce.
/// \param tag   Reduction-strategy tag (defaulted to `Pairwise{}`).
/// \returns The discrete integral as a `double`.
/// \since 0.15.0
template <class E, class Tag = Pairwise,
          std::enable_if_t<tensor::expr::is_expr_v<E> &&
                               std::is_same_v<typename std::decay_t<E>::value_type, double>,
                           int> = 0>
inline double integrate(const E& expr, Tag tag = {}) {
  static_assert(std::is_same_v<Tag, Pairwise> || std::is_same_v<Tag, Kahan>,
                "func::integrate(expr): Tag must be func::Pairwise or func::Kahan");
  return integrate(tensor::expr::materialize(expr), tag);
}

}  // namespace gridcalc::func
