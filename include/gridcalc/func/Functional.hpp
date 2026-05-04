/// \file
/// \brief End-to-end functional evaluation:
///        $F[\psi] = \int f(\psi, \nabla\psi, \nabla^2\psi, \nabla^4\psi)\, dV$.
///
/// `func::evaluate` accepts a user-supplied callable `f` and a scalar
/// `Field<double> ψ`, then materializes whatever derivatives `f` actually
/// needs (none; just `∇ψ`; both `∇ψ` and `∇²ψ`; or all four including
/// the biharmonic `∇⁴ψ`), evaluates `f` pointwise into a temporary
/// `Field<double>`, and reduces with `func::integrate`.
///
/// The callable's arity is detected at compile time via
/// `std::is_invocable_r_v`; only the derivatives the callable actually
/// consumes are materialized.
/// \since 0.4.0 (function); 0.11.0 (4-arg arity with $\nabla^4\psi$).

#pragma once

#include <type_traits>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/diff/Biharmonic.hpp>
#include <gridcalc/diff/Gradient.hpp>
#include <gridcalc/diff/Laplacian.hpp>
#include <gridcalc/func/Integrate.hpp>

namespace gridcalc::func {

namespace detail {

/// \brief Helper to defer a `static_assert(false, ...)` to template
///        instantiation, so the assert only fires when the unsupported
///        branch is actually reached.
/// \since 0.4.0
template <class>
struct deferred_false : std::false_type {};

}  // namespace detail

/// \brief Returns $F[\psi] = \int f(\psi, \nabla\psi, \nabla^2\psi, \nabla^4\psi)\, dV$.
///
/// The callable `f` may have any of four signatures (return type must be
/// convertible to `double`):
/// - `f(double psi)` — only `ψ` is needed; no derivatives are computed.
/// - `f(double psi, const core::Vec3d& grad_psi)` — `ψ` and `∇ψ` only.
/// - `f(double psi, const core::Vec3d& grad_psi, double lap_psi)` — `ψ`,
///   `∇ψ`, and the Laplacian `∇²ψ`.
/// - `f(double psi, const core::Vec3d& grad_psi, double lap_psi, double biharm_psi)`
///   — full quadruplet, with `biharm_psi` the contracted scalar
///   $\nabla^4\psi$ as returned by `diff::biharmonic` (Phase 8). Phase 11
///   added this arity for densities of up to 4th-order gradient terms
///   (e.g., the phase-field-crystal $(q_0^2+\nabla^2)^2$ kernel); rank-3
///   and rank-4 tensor-valued derivative arguments are deferred to
///   Phase 13/14.
///
/// The first matching arity is selected at compile time via
/// `std::is_invocable_r_v`; if multiple match (e.g., a callable with
/// default arguments), the longest-arity branch wins. Callables that match
/// none of the four trigger a `static_assert` at the call site.
///
/// Internally: required derivatives are materialized eagerly via
/// `diff::gradient`, `diff::laplacian`, and/or `diff::biharmonic`, the
/// callable is evaluated pointwise into a `Field<double>` integrand, and
/// the result is reduced by `func::integrate(integrand, tag)`. Phase 14
/// may replace this with expression-template / fused-loop variants once
/// tensor-valued field arguments arrive.
///
/// \tparam F     Callable type; deduced from the argument.
/// \tparam Tag   Reduction-strategy tag — must be `Pairwise` or `Kahan`.
/// \param  psi   Input scalar field.
/// \param  f     Functional integrand.
/// \param  tag   Reduction-strategy tag forwarded to `func::integrate`
///               (defaulted to `Pairwise{}`).
/// \returns The discrete functional value as a `double`.
/// \since 0.4.0 (function); 0.11.0 (4-arg arity with $\nabla^4\psi$).
template <class F, class Tag = Pairwise>
inline double evaluate(const core::Field<double>& psi, F&& f, Tag tag = {}) {
  static_assert(std::is_same_v<Tag, Pairwise> || std::is_same_v<Tag, Kahan>,
                "func::evaluate: Tag must be func::Pairwise or func::Kahan");

  const auto& grid = psi.getGrid();
  const int Nx = grid.getNx();
  const int Ny = grid.getNy();
  const int Nz = grid.getNz();

  core::Field<double> integrand(grid);

  if constexpr (std::is_invocable_r_v<double, F&, double, const core::Vec3d&, double, double>) {
    const auto grad = diff::gradient(psi);
    const auto lap = diff::laplacian(psi);
    const auto biharm = diff::biharmonic(psi);
    for (int k = 0; k < Nz; ++k) {
      for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
          integrand(i, j, k) = static_cast<double>(
              f(psi(i, j, k), grad(i, j, k), lap(i, j, k), biharm(i, j, k)));
        }
      }
    }
  } else if constexpr (std::is_invocable_r_v<double, F&, double, const core::Vec3d&, double>) {
    const auto grad = diff::gradient(psi);
    const auto lap = diff::laplacian(psi);
    for (int k = 0; k < Nz; ++k) {
      for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
          integrand(i, j, k) = static_cast<double>(f(psi(i, j, k), grad(i, j, k), lap(i, j, k)));
        }
      }
    }
  } else if constexpr (std::is_invocable_r_v<double, F&, double, const core::Vec3d&>) {
    const auto grad = diff::gradient(psi);
    for (int k = 0; k < Nz; ++k) {
      for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
          integrand(i, j, k) = static_cast<double>(f(psi(i, j, k), grad(i, j, k)));
        }
      }
    }
  } else if constexpr (std::is_invocable_r_v<double, F&, double>) {
    for (int k = 0; k < Nz; ++k) {
      for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
          integrand(i, j, k) = static_cast<double>(f(psi(i, j, k)));
        }
      }
    }
  } else {
    static_assert(detail::deferred_false<F>::value,
                  "func::evaluate: f must be invocable as f(double), "
                  "f(double, const Vec3d&), f(double, const Vec3d&, double), "
                  "or f(double, const Vec3d&, double, double), "
                  "with a return type convertible to double");
  }

  return integrate(integrand, tag);
}

}  // namespace gridcalc::func
