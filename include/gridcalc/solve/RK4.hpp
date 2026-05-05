/// \file
/// \brief Classic 4-stage Runge–Kutta time integrator for `Field<double>` state.
///
/// Defines the integrator-tag struct `gridcalc::solve::RK4` and the
/// corresponding `solve::integrate(..., RK4)` overload. Higher-order
/// time accuracy than `ExplicitEuler` at the cost of four RHS evaluations
/// per step; also has a larger heat-equation stability bound
/// (`diffusionCFLLimit = 0.6963` vs `0.5`).
/// \since 0.6.0

#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <gridcalc/core/Field.hpp>
#include <gridcalc/solve/detail/RhsGridCheck.hpp>

namespace gridcalc::solve {

/// \brief Tag selecting classic 4-stage Runge–Kutta time stepping.
///
/// Carries the heat-equation CFL bound used by `solve::diffuse`:
/// $D \cdot dt \cdot \sum_a (1/h_a^2) \le 0.6963$.
/// Derived from RK4's stability region intersecting the negative real
/// axis at $\approx -2.7853$ (Phase 6 Developer Note).
/// \since 0.6.0
struct RK4 {
  /// \brief Stability-bound coefficient for the RK4-discretized heat
  ///        equation. About 1.39× looser than `ExplicitEuler::diffusionCFLLimit`.
  /// \since 0.6.0
  static constexpr double diffusionCFLLimit = 0.6963;
};

namespace detail {

/// \brief Returns `a + scale * b` element-wise as a fresh Field.
///
/// Used to build RK4's intermediate stage inputs. `a` and `b` must
/// carry the **same `Grid`** — bit-exact equality on `(Nx, Ny, Nz)`
/// and per-axis cell sizes. The flat-index loop over `pa[n] + scale * pb[n]`
/// would otherwise silently mix cells from misaligned grids (or read
/// past the end of the smaller field's storage).
/// \param a      Base field; defines the output grid.
/// \param scale  Scalar multiplier on `b`.
/// \param b      Field to add; must share `a`'s grid.
/// \throws std::invalid_argument if `a.getGrid() != b.getGrid()`.
/// \since 0.6.0 (function); 0.14.2 (Grid-equality precondition).
inline core::Field<double> axpyFresh(const core::Field<double>& a,
                                     double scale,
                                     const core::Field<double>& b) {
  requireSameGrid(a, b, "solve::detail::axpyFresh");
  core::Field<double> out(a.getGrid());
  const std::size_t N = out.getSize();
  const double* pa = a.data();
  const double* pb = b.data();
  double* po = out.data();
  for (std::size_t n = 0; n < N; ++n) {
    po[n] = pa[n] + scale * pb[n];
  }
  return out;
}

}  // namespace detail

/// \brief Advances `psi` by `n_steps` classic-RK4 steps of size `dt`,
///        with $\partial_t\psi = \mathrm{rhs}(\psi)$.
///
/// Per step (the Butcher tableau weights `(1/6, 1/3, 1/3, 1/6)`):
/// \verbatim
///   k1 = rhs(psi)
///   k2 = rhs(psi + (dt/2) * k1)
///   k3 = rhs(psi + (dt/2) * k2)
///   k4 = rhs(psi +  dt    * k3)
///   psi += (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
/// \endverbatim
///
/// Each `psi + a * k_i` materializes one intermediate `Field<double>`;
/// each `rhs(...)` allocates one. Per-step allocation totals ~6 fields.
/// Phase 20 is expected to switch to a persistent scratch pool.
///
/// Like the `ExplicitEuler` overload, this does **not** validate
/// stability; higher-level drivers (`solve::diffuse`) consult
/// `RK4::diffusionCFLLimit` for the heat-equation bound.
/// \tparam Rhs    RHS callable type; deduced from the argument.
/// \param psi     The state field; mutated in place.
/// \param rhs     The right-hand side $\partial_t\psi = \mathrm{rhs}(\psi)$.
/// \param dt      Time step. Must be `>= 0`.
/// \param n_steps Number of steps to advance. Must be `>= 0`.
/// \param tag     Empty `RK4{}` tag for dispatch.
/// \throws std::invalid_argument if `dt < 0`, `n_steps < 0`, or any
///         per-stage `rhs(...)` call returns a `Field` allocated on a
///         different `Grid` than `psi` (the integrator's flat-index
///         accumulation requires `psi`-grid-preserving RHS values; see
///         `solve::detail::requireSameGrid`).
/// \since 0.6.0 (function); 0.14.2 (RHS-grid precondition).
template <class Rhs>
inline void integrate(core::Field<double>& psi,
                      Rhs&& rhs,
                      double dt,
                      int n_steps,
                      RK4 tag) {
  (void)tag;
  static_assert(
      std::is_invocable_r_v<core::Field<double>, Rhs&, const core::Field<double>&>,
      "solve::integrate(... RK4): rhs must be invocable as "
      "Field<double> rhs(const Field<double>&)");

  if (n_steps < 0) {
    throw std::invalid_argument(
        "solve::integrate(... RK4): n_steps must be >= 0; got " +
        std::to_string(n_steps));
  }
  if (dt < 0.0) {
    throw std::invalid_argument(
        "solve::integrate(... RK4): dt must be >= 0; got " +
        std::to_string(dt));
  }

  const std::size_t N = psi.getSize();
  const double half_dt = 0.5 * dt;
  const double sixth_dt = dt / 6.0;
  const double third_dt = dt / 3.0;

  for (int step = 0; step < n_steps; ++step) {
    // k1 at psi
    core::Field<double> k1 = rhs(static_cast<const core::Field<double>&>(psi));
    detail::requireSameGrid(psi, k1, "solve::integrate(... RK4): k1");
    // k2 at psi + (dt/2) k1
    core::Field<double> y2 = detail::axpyFresh(psi, half_dt, k1);
    core::Field<double> k2 = rhs(static_cast<const core::Field<double>&>(y2));
    detail::requireSameGrid(psi, k2, "solve::integrate(... RK4): k2");
    // k3 at psi + (dt/2) k2
    core::Field<double> y3 = detail::axpyFresh(psi, half_dt, k2);
    core::Field<double> k3 = rhs(static_cast<const core::Field<double>&>(y3));
    detail::requireSameGrid(psi, k3, "solve::integrate(... RK4): k3");
    // k4 at psi + dt k3
    core::Field<double> y4 = detail::axpyFresh(psi, dt, k3);
    core::Field<double> k4 = rhs(static_cast<const core::Field<double>&>(y4));
    detail::requireSameGrid(psi, k4, "solve::integrate(... RK4): k4");

    // psi += (dt/6) k1 + (dt/3) k2 + (dt/3) k3 + (dt/6) k4
    double* p = psi.data();
    const double* p1 = k1.data();
    const double* p2 = k2.data();
    const double* p3 = k3.data();
    const double* p4 = k4.data();
    for (std::size_t n = 0; n < N; ++n) {
      p[n] += sixth_dt * (p1[n] + p4[n]) + third_dt * (p2[n] + p3[n]);
    }
  }
}

}  // namespace gridcalc::solve
