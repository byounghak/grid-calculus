/// \file
/// \brief Rank-4 partial derivatives (15 unique partials in 3D).
///
/// Single-axis (3): `d4dx4`, `d4dy4`, `d4dz4`.
/// Two distinct axes, multiplicity (3, 1) (6): `d4dx3dy`, `d4dxdy3`,
/// `d4dx3dz`, `d4dxdz3`, `d4dy3dz`, `d4dydz3`.
/// Two distinct axes, multiplicity (2, 2) (3): `d4dx2dy2`, `d4dx2dz2`,
/// `d4dy2dz2`.
/// Three distinct axes, multiplicity (2, 1, 1) (3): `d4dx2dydz`,
/// `d4dxdy2dz`, `d4dxdydz2`.
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

/// \brief Returns $\partial^4 \psi / \partial x^4$ as a fresh same-shape field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^4 \psi / \partial x^4$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d4dx4(const core::Field<double>& psi) {
  return detail::mixedDerivative<4, 0, 0, Order>(psi);
}

/// \brief Returns $\partial^4 \psi / \partial y^4$ as a fresh same-shape field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^4 \psi / \partial y^4$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d4dy4(const core::Field<double>& psi) {
  return detail::mixedDerivative<0, 4, 0, Order>(psi);
}

/// \brief Returns $\partial^4 \psi / \partial z^4$ as a fresh same-shape field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^4 \psi / \partial z^4$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d4dz4(const core::Field<double>& psi) {
  return detail::mixedDerivative<0, 0, 4, Order>(psi);
}

/// \brief Returns $\partial^4 \psi / \partial x^3 \partial y$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^4 \psi / \partial x^3 \partial y$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d4dx3dy(const core::Field<double>& psi) {
  return detail::mixedDerivative<3, 1, 0, Order>(psi);
}

/// \brief Returns $\partial^4 \psi / \partial x \partial y^3$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^4 \psi / \partial x \partial y^3$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d4dxdy3(const core::Field<double>& psi) {
  return detail::mixedDerivative<1, 3, 0, Order>(psi);
}

/// \brief Returns $\partial^4 \psi / \partial x^3 \partial z$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^4 \psi / \partial x^3 \partial z$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d4dx3dz(const core::Field<double>& psi) {
  return detail::mixedDerivative<3, 0, 1, Order>(psi);
}

/// \brief Returns $\partial^4 \psi / \partial x \partial z^3$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^4 \psi / \partial x \partial z^3$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d4dxdz3(const core::Field<double>& psi) {
  return detail::mixedDerivative<1, 0, 3, Order>(psi);
}

/// \brief Returns $\partial^4 \psi / \partial y^3 \partial z$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^4 \psi / \partial y^3 \partial z$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d4dy3dz(const core::Field<double>& psi) {
  return detail::mixedDerivative<0, 3, 1, Order>(psi);
}

/// \brief Returns $\partial^4 \psi / \partial y \partial z^3$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^4 \psi / \partial y \partial z^3$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d4dydz3(const core::Field<double>& psi) {
  return detail::mixedDerivative<0, 1, 3, Order>(psi);
}

/// \brief Returns $\partial^4 \psi / \partial x^2 \partial y^2$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^4 \psi / \partial x^2 \partial y^2$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d4dx2dy2(const core::Field<double>& psi) {
  return detail::mixedDerivative<2, 2, 0, Order>(psi);
}

/// \brief Returns $\partial^4 \psi / \partial x^2 \partial z^2$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^4 \psi / \partial x^2 \partial z^2$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d4dx2dz2(const core::Field<double>& psi) {
  return detail::mixedDerivative<2, 0, 2, Order>(psi);
}

/// \brief Returns $\partial^4 \psi / \partial y^2 \partial z^2$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^4 \psi / \partial y^2 \partial z^2$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d4dy2dz2(const core::Field<double>& psi) {
  return detail::mixedDerivative<0, 2, 2, Order>(psi);
}

/// \brief Returns $\partial^4 \psi / \partial x^2 \partial y \partial z$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^4 \psi / \partial x^2 \partial y \partial z$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d4dx2dydz(const core::Field<double>& psi) {
  return detail::mixedDerivative<2, 1, 1, Order>(psi);
}

/// \brief Returns $\partial^4 \psi / \partial x \partial y^2 \partial z$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^4 \psi / \partial x \partial y^2 \partial z$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d4dxdy2dz(const core::Field<double>& psi) {
  return detail::mixedDerivative<1, 2, 1, Order>(psi);
}

/// \brief Returns $\partial^4 \psi / \partial x \partial y \partial z^2$ as a fresh field.
/// \tparam Order  Accuracy order (`2` default, `4` available).
/// \param psi     Input scalar field.
/// \returns A new `Field<double>` holding $\partial^4 \psi / \partial x \partial y \partial z^2$.
/// \since 0.8.0
template <int Order = 2>
inline core::Field<double> d4dxdydz2(const core::Field<double>& psi) {
  return detail::mixedDerivative<1, 1, 2, Order>(psi);
}

}  // namespace gridcalc::diff
