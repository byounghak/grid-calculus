/// \file
/// \brief Rank-2 partial derivatives: $\partial^2/\partial x_i \partial x_j$
///        for every ordered pair (six functions).
///
/// Single-axis: `dxx`, `dyy`, `dzz`. Two-axis (mixed): `dxy`, `dxz`, `dyz`.
/// Naming convention: `d<rank>d<axis-with-exponent>...` with lexicographic
/// axis order and exponent `1` omitted. All functions are templated on
/// `<int Order = 2>` and forward to
/// `gridcalc::diff::detail::mixedDerivative<Nx, Ny, Nz, Order>`.
/// \since 0.8.0

#pragma once

#include <gridcalc/core/Field.hpp>
#include <gridcalc/diff/detail/MultiIndexDerivative.hpp>

namespace gridcalc::diff {

/// \brief Returns $\partial^2 \psi / \partial x^2$ as a fresh same-shape field.
/// \tparam Order  Accuracy order (`2` default, `4` available). Unsupported
///                orders trip a compile error inside the underlying stencil
///                table.
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^2 \psi / \partial x^2$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> dxx(const core::Field<double>& psi) {
  return detail::mixedDerivative<2, 0, 0, Order>(psi);
}

/// \brief Returns $\partial^2 \psi / \partial y^2$ as a fresh same-shape field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^2 \psi / \partial y^2$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> dyy(const core::Field<double>& psi) {
  return detail::mixedDerivative<0, 2, 0, Order>(psi);
}

/// \brief Returns $\partial^2 \psi / \partial z^2$ as a fresh same-shape field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^2 \psi / \partial z^2$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> dzz(const core::Field<double>& psi) {
  return detail::mixedDerivative<0, 0, 2, Order>(psi);
}

/// \brief Returns $\partial^2 \psi / \partial x \partial y$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^2 \psi / \partial x \partial y$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> dxy(const core::Field<double>& psi) {
  return detail::mixedDerivative<1, 1, 0, Order>(psi);
}

/// \brief Returns $\partial^2 \psi / \partial x \partial z$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^2 \psi / \partial x \partial z$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> dxz(const core::Field<double>& psi) {
  return detail::mixedDerivative<1, 0, 1, Order>(psi);
}

/// \brief Returns $\partial^2 \psi / \partial y \partial z$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^2 \psi / \partial y \partial z$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> dyz(const core::Field<double>& psi) {
  return detail::mixedDerivative<0, 1, 1, Order>(psi);
}

}  // namespace gridcalc::diff
