/// \file
/// \brief Convenience driver for the heat equation
///        $\partial_t \psi = D \nabla^2 \psi$ via explicit Euler.
///
/// Forwards to `solve::explicitEuler` with `diff::laplacian` as the RHS
/// and `D * dt` as the effective step. Performs a von Neumann CFL check
/// before integration begins and throws if the step is unstable.
/// \since 0.5.0

#pragma once

#include <stdexcept>
#include <string>

#include <gridcalc/core/Field.hpp>
#include <gridcalc/diff/Laplacian.hpp>
#include <gridcalc/solve/ExplicitEuler.hpp>

namespace gridcalc::solve {

namespace detail {

/// \brief Throws `std::invalid_argument` if the explicit-diffusion CFL
///        bound is violated, otherwise returns silently.
///
/// The von Neumann stability criterion for forward-Euler discretization
/// of $\partial_t\psi = D\nabla^2\psi$ on a Cartesian grid with per-axis
/// cell sizes $(h_x, h_y, h_z)$ is
/// $D \cdot dt \cdot \left(1/h_x^2 + 1/h_y^2 + 1/h_z^2\right) \le 1/2$.
/// For isotropic spacing this reduces to $D\,dt/h^2 \le 1/(2d)$, the
/// shorthand quoted in the roadmap.
/// \param grid  The sampling grid; provides per-axis cell sizes.
/// \param D     Diffusivity. Must be `>= 0`.
/// \param dt    Time step.
/// \throws std::invalid_argument when the bound is violated, with a
///         message that names the offending values.
/// \since 0.5.0
inline void checkExplicitDiffusionCFL(const core::Grid& grid, double D, double dt) {
  const auto& h = grid.getCellSize();
  const double inv_sum = 1.0 / (h(0) * h(0)) + 1.0 / (h(1) * h(1)) + 1.0 / (h(2) * h(2));
  const double cfl = D * dt * inv_sum;
  if (cfl > 0.5) {
    throw std::invalid_argument(
        std::string("solve::diffuse: explicit-diffusion CFL violated: ") +
        "D=" + std::to_string(D) + ", dt=" + std::to_string(dt) +
        ", h=(" + std::to_string(h(0)) + ", " + std::to_string(h(1)) + ", " +
        std::to_string(h(2)) + "); D*dt*(1/hx^2 + 1/hy^2 + 1/hz^2) = " +
        std::to_string(cfl) + " > 0.5");
  }
}

}  // namespace detail

/// \brief Advances `psi` by `n_steps` forward-Euler steps of the heat
///        equation $\partial_t\psi = D\nabla^2\psi$ at time-step `dt`.
///
/// Internally calls
/// `explicitEuler(psi, &diff::laplacian, D * dt, n_steps)` — the `D`
/// factor is folded into the explicit-Euler step coefficient so each
/// step computes $\psi^{n+1} = \psi^n + (D \cdot dt)\,\nabla^2 \psi^n$,
/// which is mathematically identical to running an unscaled Laplacian
/// RHS at a scaled time step. The user-visible `dt` and `n_steps` retain
/// their physical meaning: total integration time is `n_steps * dt`.
///
/// Stability is enforced before the loop runs: the von Neumann bound
/// $D \cdot dt \cdot \sum_a (1/h_a^2) \le 1/2$ must hold or the function
/// throws.
/// \param psi      The state field; mutated in place.
/// \param D        Diffusivity coefficient. Must be `>= 0`.
/// \param dt       Time step. Must be `>= 0`.
/// \param n_steps  Number of explicit-Euler steps to take. Must be `>= 0`.
/// \throws std::invalid_argument on negative `D`, `dt`, or `n_steps`,
///         or on CFL violation.
/// \since 0.5.0
inline void diffuse(core::Field<double>& psi, double D, double dt, int n_steps) {
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
  detail::checkExplicitDiffusionCFL(psi.getGrid(), D, dt);
  explicitEuler(
      psi,
      [](const core::Field<double>& f) { return diff::laplacian(f); },
      D * dt,
      n_steps);
}

}  // namespace gridcalc::solve
