/// \file
/// \brief Forward (explicit) Euler time integrator for `Field<double>` state.
///
/// Defines the integrator-tag struct `gridcalc::solve::ExplicitEuler` and
/// the corresponding `solve::integrate(..., ExplicitEuler)` overload that
/// advances $\psi^{n+1} = \psi^n + dt \cdot \mathrm{rhs}(\psi^n)$ for a
/// user-supplied right-hand-side callable.
///
/// The Phase-5 free function `solve::explicitEuler(psi, rhs, dt, n_steps)`
/// is retained as a thin wrapper that forwards to
/// `solve::integrate(psi, rhs, dt, n_steps, ExplicitEuler{})`.
///
/// Stability is **not** checked here — the integrator cannot know the
/// relevant CFL bound for an arbitrary RHS. Higher-level drivers like
/// `solve::diffuse` use `ExplicitEuler::diffusionCFLLimit` to enforce the
/// heat-equation CFL.
/// \since 0.5.0

#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <gridcalc/core/Field.hpp>

namespace gridcalc::solve {

/// \brief Tag selecting forward (explicit) Euler time stepping.
///
/// Carries the heat-equation CFL bound used by `solve::diffuse`:
/// $D \cdot dt \cdot \sum_a (1/h_a^2) \le 0.5$.
/// Derived from von-Neumann analysis (Phase 5 Developer Note).
/// \since 0.6.0
struct ExplicitEuler {
  /// \brief Stability-bound coefficient for the explicit-Euler-discretized
  ///        heat equation. \since 0.6.0
  static constexpr double diffusionCFLLimit = 0.5;
};

namespace detail {

/// \brief Helper to defer a `static_assert(false, ...)` to template
///        instantiation, matching Phase 4's `detail::deferred_false`
///        trick under C++17.
/// \since 0.5.0
template <class>
struct deferred_false : std::false_type {};

}  // namespace detail

/// \brief Advances `psi` by `n_steps` forward-Euler steps of size `dt`,
///        with $\partial_t\psi = \mathrm{rhs}(\psi)$.
///
/// At each step:
/// - `tmp = rhs(psi)` is evaluated (the RHS callable returns a fresh
///   `Field<double>` of the same shape).
/// - `psi[n] += dt * tmp[n]` element-wise via raw pointer access.
///
/// The integrator does **not** validate stability; the caller (or a
/// higher-level driver such as `solve::diffuse`) is responsible for any
/// CFL-style argument that depends on the specific RHS.
/// \tparam Rhs    RHS callable type; deduced from the argument.
/// \param psi     The state field; mutated in place.
/// \param rhs     The right-hand side $\partial_t\psi = \mathrm{rhs}(\psi)$.
///                Must be invocable as
///                `Field<double> rhs(const Field<double>&)`.
/// \param dt      Time step. Must be `>= 0`.
/// \param n_steps Number of steps to advance. Must be `>= 0`.
/// \param tag     Empty `ExplicitEuler{}` tag for dispatch.
/// \throws std::invalid_argument if `dt < 0` or `n_steps < 0`.
/// \since 0.6.0
template <class Rhs>
inline void integrate(core::Field<double>& psi,
                      Rhs&& rhs,
                      double dt,
                      int n_steps,
                      ExplicitEuler tag) {
  (void)tag;
  static_assert(
      std::is_invocable_r_v<core::Field<double>, Rhs&, const core::Field<double>&>,
      "solve::integrate(... ExplicitEuler): rhs must be invocable as "
      "Field<double> rhs(const Field<double>&)");

  if (n_steps < 0) {
    throw std::invalid_argument(
        "solve::integrate(... ExplicitEuler): n_steps must be >= 0; got " +
        std::to_string(n_steps));
  }
  if (dt < 0.0) {
    throw std::invalid_argument(
        "solve::integrate(... ExplicitEuler): dt must be >= 0; got " +
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

/// \brief Phase-5 thin wrapper. Equivalent to
///        `solve::integrate(psi, rhs, dt, n_steps, ExplicitEuler{})`.
///
/// Retained for backward compatibility — Phase 5 callers continue to
/// compile and behave identically. Prefer the generic
/// `solve::integrate(...)` form for new code, especially when callers
/// may want to swap integrators.
/// \tparam Rhs   RHS callable type.
/// \param psi    See `integrate(...)`.
/// \param rhs    See `integrate(...)`.
/// \param dt     See `integrate(...)`.
/// \param n_steps See `integrate(...)`.
/// \throws std::invalid_argument as for `integrate(... ExplicitEuler{})`.
/// \since 0.5.0
template <class Rhs>
inline void explicitEuler(core::Field<double>& psi,
                          Rhs&& rhs,
                          double dt,
                          int n_steps) {
  integrate(psi, std::forward<Rhs>(rhs), dt, n_steps, ExplicitEuler{});
}

}  // namespace gridcalc::solve
