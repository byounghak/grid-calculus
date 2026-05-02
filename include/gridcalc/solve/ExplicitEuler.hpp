/// \file
/// \brief Generic forward-Euler time integrator for `Field<double>` state.
///
/// Advances $\psi^{n+1} = \psi^n + dt \cdot \mathrm{rhs}(\psi^n)$ for a
/// user-supplied right-hand-side callable. The integrator does not check
/// stability (it cannot — that depends on the RHS); higher-level drivers
/// like `solve::diffuse` are responsible for any RHS-specific stability
/// argument and the corresponding CFL check.
/// \since 0.5.0

#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <gridcalc/core/Field.hpp>

namespace gridcalc::solve {

namespace detail {

/// \brief Helper to defer a `static_assert(false, ...)` to template
///        instantiation, matching Phase 4's `detail::deferred_false`
///        trick under C++17.
/// \since 0.5.0
template <class>
struct deferred_false : std::false_type {};

}  // namespace detail

/// \brief Advances `psi` by `n_steps` forward-Euler steps of size `dt`,
///        with $\dot\psi = \mathrm{rhs}(\psi)$.
///
/// At each step:
/// - `tmp = rhs(psi)` is evaluated (the RHS callable returns a fresh
///   `Field<double>` of the same shape).
/// - `psi[n] += dt * tmp[n]` element-wise via raw pointer access.
///
/// The integrator does **not** validate stability; the caller (or a
/// higher-level driver such as `solve::diffuse`) is responsible for
/// any CFL-style argument that depends on the specific RHS.
/// \tparam Rhs    RHS callable type; deduced from the argument.
/// \param psi     The state field; mutated in place. Updated to the
///                value at time `t + n_steps * dt`.
/// \param rhs     The right-hand side $\dot\psi = \mathrm{rhs}(\psi)$.
///                Must be invocable as
///                `Field<double> rhs(const Field<double>&)` (return type
///                must be convertible to `Field<double>`).
/// \param dt      Time step. Must be `>= 0`.
/// \param n_steps Number of steps to advance. Must be `>= 0`.
/// \throws std::invalid_argument if `dt < 0` or `n_steps < 0`.
/// \since 0.5.0
template <class Rhs>
inline void explicitEuler(core::Field<double>& psi, Rhs&& rhs, double dt, int n_steps) {
  static_assert(
      std::is_invocable_r_v<core::Field<double>, Rhs&, const core::Field<double>&>,
      "solve::explicitEuler: rhs must be invocable as "
      "Field<double> rhs(const Field<double>&)");

  if (n_steps < 0) {
    throw std::invalid_argument("solve::explicitEuler: n_steps must be >= 0; got " +
                                std::to_string(n_steps));
  }
  if (dt < 0.0) {
    throw std::invalid_argument("solve::explicitEuler: dt must be >= 0; got " +
                                std::to_string(dt));
  }

  const std::size_t N = psi.getSize();
  for (int step = 0; step < n_steps; ++step) {
    core::Field<double> tmp = rhs(static_cast<const core::Field<double>&>(psi));
    double* p = psi.data();
    const double* r = tmp.data();
    for (std::size_t n = 0; n < N; ++n) {
      p[n] += dt * r[n];
    }
  }
}

}  // namespace gridcalc::solve
