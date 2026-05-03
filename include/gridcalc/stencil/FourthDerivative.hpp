/// \file
/// \brief Central-difference stencil coefficients for the **fourth**
///        derivative in 1D.
///
/// Sibling table to `gridcalc::stencil::FirstDerivative<Order>` (first
/// derivative), `gridcalc::stencil::Coefficients<Order>` (second derivative),
/// and `gridcalc::stencil::ThirdDerivative<Order>` (third derivative).
/// Phase 8 ships `Order = 2` and `Order = 4` specializations. The primary
/// template is intentionally undefined so that an unsupported `Order`
/// produces a compile error at the call site rather than a silent
/// "all zero" stencil.
/// \since 0.8.0

#pragma once

#include <array>

namespace gridcalc::stencil {

/// \brief Central-difference coefficients for the **fourth derivative** at
///        accuracy order `Order`.
///
/// Specializations expose:
/// - `static constexpr int radius` — half-width of the stencil; the stencil
///   reaches indices in `[-radius, +radius]`.
/// - `static constexpr std::array<int, 2*radius+1> offsets`
/// - `static constexpr std::array<double, 2*radius+1> weights`
///
/// Weights are **unscaled**: the consumer multiplies the stencil sum by
/// `1 / h^4` to obtain $\partial^4 / \partial x^4$. This separation lets the
/// same coefficient table drive operators with different cell sizes per axis
/// (`Field` carries `Vec3d {hx, hy, hz}`).
///
/// The primary template is intentionally undefined; instantiating an
/// unsupported `Order` produces a compile error pointing at the call site.
/// \tparam Order  Accuracy order of the central-difference stencil (must be
///                a supported even integer; 2 and 4 at Phase 8).
/// \since 0.8.0
template <int Order>
struct FourthDerivative;

/// \brief 2nd-order central-difference coefficients for $\partial^4/\partial x^4$.
///
/// Stencil: $\partial^4 f / \partial x^4 \approx
///   (f_{-2} - 4 f_{-1} + 6 f_0 - 4 f_{+1} + f_{+2}) / h^4$.
/// Truncation error is $O(h^2)$. The derivation by Taylor matching is in
/// the Phase 8 Developer Note.
/// \since 0.8.0
template <>
struct FourthDerivative<2> {
  /// \brief Stencil half-width; the stencil reaches indices in `[-2, +2]`.
  /// \since 0.8.0
  static constexpr int radius = 2;
  /// \brief Stencil offsets relative to the evaluation point.
  /// \since 0.8.0
  static constexpr std::array<int, 5> offsets = {-2, -1, 0, 1, 2};
  /// \brief Stencil weights, unscaled (multiply by `1 / h^4` at the consumer).
  /// \since 0.8.0
  static constexpr std::array<double, 5> weights = {1.0, -4.0, 6.0, -4.0, 1.0};
};

/// \brief 4th-order central-difference coefficients for $\partial^4/\partial x^4$.
///
/// Stencil: $\partial^4 f / \partial x^4 \approx
///   (-\tfrac{1}{6} f_{-3} + 2 f_{-2} - \tfrac{13}{2} f_{-1}
///    + \tfrac{28}{3} f_0
///    - \tfrac{13}{2} f_{+1} + 2 f_{+2} - \tfrac{1}{6} f_{+3}) / h^4$.
/// Truncation error is $O(h^4)$. The derivation by Taylor matching is in
/// the Phase 8 Developer Note.
/// \since 0.8.0
template <>
struct FourthDerivative<4> {
  /// \brief Stencil half-width; the stencil reaches indices in `[-3, +3]`.
  /// \since 0.8.0
  static constexpr int radius = 3;
  /// \brief Stencil offsets relative to the evaluation point.
  /// \since 0.8.0
  static constexpr std::array<int, 7> offsets = {-3, -2, -1, 0, 1, 2, 3};
  /// \brief Stencil weights, unscaled (multiply by `1 / h^4` at the consumer).
  /// \since 0.8.0
  static constexpr std::array<double, 7> weights = {
      -1.0 / 6.0, 2.0, -13.0 / 2.0, 28.0 / 3.0, -13.0 / 2.0, 2.0, -1.0 / 6.0};
};

}  // namespace gridcalc::stencil
