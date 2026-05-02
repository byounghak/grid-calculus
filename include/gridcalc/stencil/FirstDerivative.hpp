/// \file
/// \brief Central-difference stencil coefficients for the **first**
///        derivative in 1D.
///
/// Sibling table to `gridcalc::stencil::Coefficients<Order>` (which holds
/// **second**-derivative coefficients). At Phase 2 only the 2nd-order
/// specialization is provided; Phase 7 adds `Order = 4` alongside
/// `Coefficients<4>`.
///
/// The two derivative-rank tables are intentionally separate templates with
/// distinct names so that Phase 1's public contract for `Coefficients<Order>`
/// (= second derivative) is preserved. A unifying redesign — for example a
/// 2-parameter `Coefficients<Derivative, Order>` — is deferred to Phase 7
/// when both tables gain `Order = 4` specializations and the trade-offs are
/// easier to evaluate.
/// \since 0.2.0

#pragma once

#include <array>

namespace gridcalc::stencil {

/// \brief Central-difference coefficients for the **first derivative** at
///        accuracy order `Order`.
///
/// Specializations expose:
/// - `static constexpr int radius` — half-width of the stencil; the stencil
///   reaches indices in `[-radius, +radius]`.
/// - `static constexpr std::array<int, 2*radius+1> offsets`
/// - `static constexpr std::array<double, 2*radius+1> weights`
///
/// Weights are **unscaled**: the consumer multiplies the stencil sum by
/// `1 / h` to obtain $\partial / \partial x$. This separation lets the same
/// coefficient table drive operators with different cell sizes per axis
/// (`Field` carries `Vec3d {hx, hy, hz}`).
///
/// The primary template is intentionally undefined; instantiating an
/// unsupported `Order` produces a compile error pointing at the call site
/// rather than a silent "all zero" stencil.
/// \tparam Order  Accuracy order of the central-difference stencil (must be
///                a supported even integer; 2 at Phase 2).
/// \since 0.2.0
template <int Order>
struct FirstDerivative;

/// \brief 2nd-order central-difference coefficients for $\partial/\partial x$.
///
/// Stencil: $\partial f / \partial x \approx (f_{+1} - f_{-1}) / (2 h)$.
/// \since 0.2.0
template <>
struct FirstDerivative<2> {
  /// \brief Stencil half-width; the stencil reaches indices in `[-1, +1]`.
  /// \since 0.2.0
  static constexpr int radius = 1;
  /// \brief Stencil offsets relative to the evaluation point.
  /// \since 0.2.0
  static constexpr std::array<int, 3> offsets = {-1, 0, 1};
  /// \brief Stencil weights, unscaled (multiply by `1 / h` at the consumer).
  /// \since 0.2.0
  static constexpr std::array<double, 3> weights = {-0.5, 0.0, 0.5};
};

}  // namespace gridcalc::stencil
