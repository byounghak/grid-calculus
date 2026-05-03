/// \file
/// \brief Central-difference stencil coefficients for the **third**
///        derivative in 1D.
///
/// Sibling table to `gridcalc::stencil::FirstDerivative<Order>` (first
/// derivative) and `gridcalc::stencil::Coefficients<Order>` (second
/// derivative). Phase 8 ships `Order = 2` and `Order = 4` specializations.
/// Like its siblings the primary template is intentionally undefined so
/// that an unsupported `Order` produces a compile error at the call site
/// rather than a silent "all zero" stencil.
/// \since 0.8.0

#pragma once

#include <array>

namespace gridcalc::stencil {

/// \brief Central-difference coefficients for the **third derivative** at
///        accuracy order `Order`.
///
/// Specializations expose:
/// - `static constexpr int radius` — half-width of the stencil; the stencil
///   reaches indices in `[-radius, +radius]`.
/// - `static constexpr std::array<int, 2*radius+1> offsets`
/// - `static constexpr std::array<double, 2*radius+1> weights`
///
/// Weights are **unscaled**: the consumer multiplies the stencil sum by
/// `1 / h^3` to obtain $\partial^3 / \partial x^3$. This separation lets the
/// same coefficient table drive operators with different cell sizes per axis
/// (`Field` carries `Vec3d {hx, hy, hz}`).
///
/// The primary template is intentionally undefined; instantiating an
/// unsupported `Order` produces a compile error pointing at the call site.
/// \tparam Order  Accuracy order of the central-difference stencil (must be
///                a supported even integer; 2 and 4 at Phase 8).
/// \since 0.8.0
template <int Order>
struct ThirdDerivative;

/// \brief 2nd-order central-difference coefficients for $\partial^3/\partial x^3$.
///
/// Stencil: $\partial^3 f / \partial x^3 \approx
///   (-\tfrac{1}{2} f_{-2} + f_{-1} - f_{+1} + \tfrac{1}{2} f_{+2}) / h^3$.
/// Truncation error is $O(h^2)$. The center coefficient is exactly zero
/// (anti-symmetric stencil for an odd derivative). The derivation by Taylor
/// matching is in the Phase 8 Developer Note.
/// \since 0.8.0
template <>
struct ThirdDerivative<2> {
  /// \brief Stencil half-width; the stencil reaches indices in `[-2, +2]`.
  /// \since 0.8.0
  static constexpr int radius = 2;
  /// \brief Stencil offsets relative to the evaluation point.
  /// \since 0.8.0
  static constexpr std::array<int, 5> offsets = {-2, -1, 0, 1, 2};
  /// \brief Stencil weights, unscaled (multiply by `1 / h^3` at the consumer).
  /// \since 0.8.0
  static constexpr std::array<double, 5> weights = {
      -1.0 / 2.0, 1.0, 0.0, -1.0, 1.0 / 2.0};
};

/// \brief 4th-order central-difference coefficients for $\partial^3/\partial x^3$.
///
/// Stencil: $\partial^3 f / \partial x^3 \approx
///   (\tfrac{1}{8} f_{-3} - f_{-2} + \tfrac{13}{8} f_{-1}
///    - \tfrac{13}{8} f_{+1} + f_{+2} - \tfrac{1}{8} f_{+3}) / h^3$.
/// Truncation error is $O(h^4)$. The center coefficient is exactly zero
/// (anti-symmetric stencil for an odd derivative). The derivation by Taylor
/// matching is in the Phase 8 Developer Note.
/// \since 0.8.0
template <>
struct ThirdDerivative<4> {
  /// \brief Stencil half-width; the stencil reaches indices in `[-3, +3]`.
  /// \since 0.8.0
  static constexpr int radius = 3;
  /// \brief Stencil offsets relative to the evaluation point.
  /// \since 0.8.0
  static constexpr std::array<int, 7> offsets = {-3, -2, -1, 0, 1, 2, 3};
  /// \brief Stencil weights, unscaled (multiply by `1 / h^3` at the consumer).
  /// \since 0.8.0
  static constexpr std::array<double, 7> weights = {
      1.0 / 8.0, -1.0, 13.0 / 8.0, 0.0, -13.0 / 8.0, 1.0, -1.0 / 8.0};
};

}  // namespace gridcalc::stencil
