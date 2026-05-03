/// \file
/// \brief Rank-3 partial derivatives (10 unique partials in 3D).
///
/// Single-axis: `d3dx3`, `d3dy3`, `d3dz3`.
/// Two distinct axes (multiplicity (2, 1)): `d3dx2dy`, `d3dxdy2`,
/// `d3dx2dz`, `d3dxdz2`, `d3dy2dz`, `d3dydz2`.
/// Three distinct axes: `d3dxdydz`.
///
/// Naming convention: `d<rank>d<axis-with-exponent>...` with lexicographic
/// axis order and exponent `1` omitted. All functions are templated on
/// `<int Order = 2>` and forward to
/// `gridcalc::diff::detail::mixedDerivative<Nx, Ny, Nz, Order>`.
/// \since 0.8.0

#pragma once

#include <gridcalc/core/Field.hpp>
#include <gridcalc/diff/detail/MultiIndexDerivative.hpp>

namespace gridcalc::diff {

/// \brief Returns $\partial^3 \psi / \partial x^3$ as a fresh same-shape field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^3 \psi / \partial x^3$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d3dx3(const core::Field<double>& psi) {
  return detail::mixedDerivative<3, 0, 0, Order>(psi);
}

/// \brief Returns $\partial^3 \psi / \partial y^3$ as a fresh same-shape field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^3 \psi / \partial y^3$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d3dy3(const core::Field<double>& psi) {
  return detail::mixedDerivative<0, 3, 0, Order>(psi);
}

/// \brief Returns $\partial^3 \psi / \partial z^3$ as a fresh same-shape field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^3 \psi / \partial z^3$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d3dz3(const core::Field<double>& psi) {
  return detail::mixedDerivative<0, 0, 3, Order>(psi);
}

/// \brief Returns $\partial^3 \psi / \partial x^2 \partial y$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^3 \psi / \partial x^2 \partial y$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d3dx2dy(const core::Field<double>& psi) {
  return detail::mixedDerivative<2, 1, 0, Order>(psi);
}

/// \brief Returns $\partial^3 \psi / \partial x \partial y^2$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^3 \psi / \partial x \partial y^2$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d3dxdy2(const core::Field<double>& psi) {
  return detail::mixedDerivative<1, 2, 0, Order>(psi);
}

/// \brief Returns $\partial^3 \psi / \partial x^2 \partial z$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^3 \psi / \partial x^2 \partial z$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d3dx2dz(const core::Field<double>& psi) {
  return detail::mixedDerivative<2, 0, 1, Order>(psi);
}

/// \brief Returns $\partial^3 \psi / \partial x \partial z^2$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^3 \psi / \partial x \partial z^2$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d3dxdz2(const core::Field<double>& psi) {
  return detail::mixedDerivative<1, 0, 2, Order>(psi);
}

/// \brief Returns $\partial^3 \psi / \partial y^2 \partial z$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^3 \psi / \partial y^2 \partial z$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d3dy2dz(const core::Field<double>& psi) {
  return detail::mixedDerivative<0, 2, 1, Order>(psi);
}

/// \brief Returns $\partial^3 \psi / \partial y \partial z^2$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^3 \psi / \partial y \partial z^2$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d3dydz2(const core::Field<double>& psi) {
  return detail::mixedDerivative<0, 1, 2, Order>(psi);
}

/// \brief Returns $\partial^3 \psi / \partial x \partial y \partial z$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^3 \psi / \partial x \partial y \partial z$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d3dxdydz(const core::Field<double>& psi) {
  return detail::mixedDerivative<1, 1, 1, Order>(psi);
}

}  // namespace gridcalc::diff
