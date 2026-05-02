/// \file
/// \brief Central-difference stencil coefficients for derivatives in 1D.
///
/// At Phase 1 only the 2nd-order specialization for the second derivative is
/// provided. Phase 7 adds `Order = 4`.
/// \since 0.1.0

#pragma once

#include <array>

namespace gridcalc::stencil {

/// \brief Central-difference coefficients for the **second derivative** at
///        accuracy order `Order`.
///
/// Specializations expose:
/// - `static constexpr int radius` — half-width of the stencil; the stencil
///   reaches indices in `[-radius, +radius]`.
/// - `static constexpr std::array<int, 2*radius+1> offsets`
/// - `static constexpr std::array<double, 2*radius+1> weights`
///
/// Weights are **unscaled**: the consumer multiplies the stencil sum by
/// `1 / h^2` to obtain $\partial^2 / \partial x^2$. This separation lets the
/// same coefficient tables drive operators with different cell sizes per axis
/// (`Field` carries `Vec3d {hx, hy, hz}`).
///
/// The primary template is intentionally undefined; instantiating an
/// unsupported `Order` produces a compile error pointing at the call site
/// rather than a silent "all zero" stencil.
/// \tparam Order  Accuracy order of the central-difference stencil (must be
///                a supported even integer; 2 at Phase 1).
/// \since 0.1.0
template <int Order>
struct Coefficients;

/// \brief 2nd-order central-difference coefficients for $\partial^2/\partial x^2$.
///
/// Stencil: $\partial^2 f / \partial x^2 \approx (f_{-1} - 2 f_0 + f_{+1}) / h^2$.
/// \since 0.1.0
template <>
struct Coefficients<2> {
  /// \brief Stencil half-width; the stencil reaches indices in `[-1, +1]`.
  /// \since 0.1.0
  static constexpr int radius = 1;
  /// \brief Stencil offsets relative to the evaluation point.
  /// \since 0.1.0
  static constexpr std::array<int, 3> offsets = {-1, 0, 1};
  /// \brief Stencil weights, unscaled (multiply by `1 / h^2` at the consumer).
  /// \since 0.1.0
  static constexpr std::array<double, 3> weights = {1.0, -2.0, 1.0};
};

}  // namespace gridcalc::stencil
