/// \file
/// \brief Convenience driver for the heat equation
///        $\partial_t \psi = D \nabla^2 \psi$ via a chosen time integrator.
///
/// `solve::diffuse(psi, D, dt, n_steps)` defaults to forward Euler
/// (Phase 5 behavior). Pass an explicit `RK4{}` tag to switch:
/// `solve::diffuse(psi, D, dt, n_steps, RK4{})`. The CFL stability
/// check reads `Tag::diffusionCFLLimit`, so each integrator gets the
/// stability bound appropriate to its discrete amplification factor.
/// \since 0.5.0

#pragma once

#include <algorithm>
#include <stdexcept>
#include <string>
#include <typeinfo>

#include <gridcalc/core/Field.hpp>
#include <gridcalc/diff/Laplacian.hpp>
#include <gridcalc/fvm/CellLaplacian.hpp>
#include <gridcalc/solve/Integrator.hpp>

namespace gridcalc::solve {

namespace detail {

/// \brief Throws `std::invalid_argument` if the explicit-diffusion CFL
///        bound is violated for the given integrator tag, otherwise
///        returns silently.
///
/// The von-Neumann criterion for forward-Euler discretization of
/// $\partial_t\psi = D\nabla^2\psi$ on a Cartesian grid is
/// $D \cdot dt \cdot \sum_a (1/h_a^2) \le 0.5$. RK4 admits a larger
/// bound (`~0.6963`) thanks to its stability region. The integrator's
/// tag struct exposes its own bound as `Tag::diffusionCFLLimit`.
/// \tparam Tag  Integrator tag exposing a `static constexpr double diffusionCFLLimit`.
/// \param grid  The sampling grid; provides per-axis cell sizes.
/// \param D     Diffusivity. Must be `>= 0`.
/// \param dt    Time step.
/// \throws std::invalid_argument when the bound is violated, with a
///         message that names the offending values and limit.
/// \since 0.6.0
template <class Tag>
inline void checkDiffusionCFL(const core::Grid& grid, double D, double dt) {
  const auto& h = grid.getCellSize();
  const double inv_sum = 1.0 / (h(0) * h(0)) + 1.0 / (h(1) * h(1)) + 1.0 / (h(2) * h(2));
  const double cfl = D * dt * inv_sum;
  if (cfl > Tag::diffusionCFLLimit) {
    throw std::invalid_argument(
        std::string("solve::diffuse: explicit-diffusion CFL violated for ") +
        typeid(Tag).name() +
        ": D=" + std::to_string(D) + ", dt=" + std::to_string(dt) +
        ", h=(" + std::to_string(h(0)) + ", " + std::to_string(h(1)) + ", " +
        std::to_string(h(2)) + "); D*dt*(1/hx^2 + 1/hy^2 + 1/hz^2) = " +
        std::to_string(cfl) + " > " + std::to_string(Tag::diffusionCFLLimit));
  }
}

}  // namespace detail

/// \brief Advances `psi` by `n_steps` steps of the heat equation
///        $\partial_t\psi = D\nabla^2\psi$ at time-step `dt`, using
///        the integrator selected by `Integrator`.
///
/// Internally calls
/// `solve::integrate(psi, &diff::laplacian, D * dt, n_steps, tag)`
/// — the `D` factor is folded into the integrator's effective step
/// (mathematically identical to running an unscaled Laplacian RHS at a
/// scaled time step). The user-visible `dt` and `n_steps` retain their
/// physical meaning.
///
/// Stability is enforced before the loop runs: the von-Neumann bound
/// `D · dt · sum_a(1/h_a^2) <= Integrator::diffusionCFLLimit` must
/// hold or the function throws.
/// \tparam Integrator  Integrator tag (`ExplicitEuler` default, or `RK4`).
/// \param psi      The state field; mutated in place.
/// \param D        Diffusivity coefficient. Must be `>= 0`.
/// \param dt       Time step. Must be `>= 0`.
/// \param n_steps  Number of integrator steps to take. Must be `>= 0`.
/// \param tag      Integrator selector (defaulted to `ExplicitEuler{}`).
/// \throws std::invalid_argument on negative `D`, `dt`, or `n_steps`,
///         on CFL violation for the selected integrator, or on a grid
///         too small for `diff::laplacian`'s stencil radius (propagated
///         from the underlying operator; see
///         `diff::detail::requireAxisExtent`).
/// \since 0.5.0
template <class Integrator = ExplicitEuler>
inline void diffuse(core::Field<double>& psi,
                    double D,
                    double dt,
                    int n_steps,
                    Integrator tag = {}) {
  if (D < 0.0) {
    throw std::invalid_argument("solve::diffuse: D must be >= 0; got " + std::to_string(D));
  }
  if (dt < 0.0) {
    throw std::invalid_argument("solve::diffuse: dt must be >= 0; got " + std::to_string(dt));
  }
  if (n_steps < 0) {
    throw std::invalid_argument("solve::diffuse: n_steps must be >= 0; got " +
                                std::to_string(n_steps));
  }
  detail::checkDiffusionCFL<Integrator>(psi.getGrid(), D, dt);
  integrate(
      psi,
      [](const core::Field<double>& f) { return diff::laplacian(f); },
      D * dt,
      n_steps,
      tag);
}

/// \brief Advances `psi` by `n_steps` steps of the **heterogeneous-D**
///        diffusion equation $\partial_t\psi = \nabla\cdot(D\,\nabla\psi)$
///        at time-step `dt`, using the integrator selected by `Integrator`.
///
/// Routes the RHS through `fvm::cellLaplacian(psi, D)` (Phase 14
/// finite-volume cell-flux discretization with face-centered
/// harmonic-mean D-averaging). The constant-D Phase 5 / Phase 6 overload
/// (above) is unchanged and stays the right choice for callers with a
/// uniform diffusivity — on a uniform Cartesian grid both paths
/// produce the same trajectory (verified by
/// `HeterogeneousD_AgreesWithConstantDOnUniformField`).
///
/// **CFL bound.** The von-Neumann criterion generalizes naturally:
/// $D_{\max}\cdot dt\cdot\sum_a (1/h_a^2)\le \texttt{Tag::diffusionCFLLimit}$,
/// where $D_{\max}$ is the spatial maximum of `D`. Stricter than the
/// constant-D bound when `D` has high-conductivity spots; equal to it
/// when `D` is uniform.
///
/// \tparam Integrator  Integrator tag (`ExplicitEuler` default, or `RK4`).
/// \param psi      The state field; mutated in place.
/// \param D        Diffusivity field on the same `Grid` as `psi`. Must
///                 be strictly positive at every cell (the harmonic-mean
///                 face-D averaging used by `fvm::cellLaplacian` is
///                 undefined otherwise).
/// \param dt       Time step. Must be `>= 0`.
/// \param n_steps  Number of integrator steps to take. Must be `>= 0`.
/// \param tag      Integrator selector (defaulted to `ExplicitEuler{}`).
/// \throws std::invalid_argument on negative `dt` or `n_steps`, on
///         non-positive `D` (when caught: a representative cell that
///         violated the contract is named in the message), on CFL
///         violation for the selected integrator, or on a grid too small
///         for `fvm::cellLaplacian`'s stencil radius (propagated from
///         the underlying operator; see
///         `diff::detail::requireAxisExtent`).
/// \since 0.14.0
template <class Integrator = ExplicitEuler>
inline void diffuse(core::Field<double>& psi,
                    const core::Field<double>& D,
                    double dt,
                    int n_steps,
                    Integrator tag = {}) {
  if (dt < 0.0) {
    throw std::invalid_argument("solve::diffuse: dt must be >= 0; got " +
                                std::to_string(dt));
  }
  if (n_steps < 0) {
    throw std::invalid_argument("solve::diffuse: n_steps must be >= 0; got " +
                                std::to_string(n_steps));
  }

  // D must be strictly positive everywhere (harmonic-mean contract).
  // Also compute D_max for the CFL bound in the same pass.
  double D_max = 0.0;
  const std::size_t n = D.getSize();
  const double* pD = D.data();
  for (std::size_t idx = 0; idx < n; ++idx) {
    if (!(pD[idx] > 0.0)) {
      throw std::invalid_argument(
          "solve::diffuse(heterogeneous-D): D must be strictly positive at "
          "every cell; first violation at flat index " +
          std::to_string(idx) + " with D=" + std::to_string(pD[idx]));
    }
    if (pD[idx] > D_max) D_max = pD[idx];
  }

  detail::checkDiffusionCFL<Integrator>(psi.getGrid(), D_max, dt);
  integrate(
      psi,
      [&D](const core::Field<double>& f) { return fvm::cellLaplacian(f, D); },
      dt,
      n_steps,
      tag);
}

}  // namespace gridcalc::solve
