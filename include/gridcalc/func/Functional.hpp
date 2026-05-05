/// \file
/// \brief End-to-end functional evaluation:
///        $F[\psi] = \int f(\psi, \nabla\psi, \nabla^2\psi, \nabla^4\psi)\, dV$
///        on scalar inputs, plus
///        $F[\mathbf{u}] = \int f(\mathbf{u}, \nabla\mathbf{u})\, dV$
///        on `Field<core::Vec3d>` inputs (Phase 15).
///
/// `func::evaluate` accepts a user-supplied callable `f` and a scalar
/// `Field<double> ψ`, then materializes whatever derivatives `f` actually
/// needs (none; just `∇ψ`; both `∇ψ` and `∇²ψ`; or all four including
/// the biharmonic `∇⁴ψ`), evaluates `f` pointwise into a temporary
/// `Field<double>`, and reduces with `func::integrate`.
///
/// Phase 15 adds two `Field<core::Vec3d>` overloads with arities
/// `f(const Vec3d& u)` and `f(const Vec3d& u, const Mat3d& grad_u)` --- the
/// surface that drives the linear-elastic-energy worked example. The
/// rank-2 gradient convention follows Phase 13: `M(i, j) = ∂_j u_i`, so
/// `trace(M) = divergence(u)`. Higher derivatives of vector / tensor
/// fields are deferred to a later phase along with the rank-3 / rank-4
/// tensor types they would consume.
///
/// The callable's arity is detected at compile time via
/// `std::is_invocable_r_v`; only the derivatives the callable actually
/// consumes are materialized.
/// \since 0.4.0 (scalar function); 0.11.0 (4-arg scalar arity with
/// $\nabla^4\psi$); 0.15.0 (`Field<core::Vec3d>` overloads).

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

/// \brief Returns $F[\mathbf{u}] = \int f(\mathbf{u}, \nabla\mathbf{u})\, dV$
///        for a `core::Vec3d`-valued field input.
///
/// The callable `f` may have either of two signatures (return type must be
/// convertible to `double`):
/// - `f(const core::Vec3d& u)` --- only `\mathbf{u}` is needed; no
///   derivative is computed.
/// - `f(const core::Vec3d& u, const core::Mat3d& grad_u)` --- `\mathbf{u}`
///   and its rank-2 gradient `M(i, j) = \partial_j u_i` (Phase 13
///   convention; `trace(M) = divergence(u)` by construction).
///
/// The first matching arity is selected at compile time via
/// `std::is_invocable_r_v`; if the callable matches both (e.g., a generic
/// lambda taking `auto`), the longest-arity branch wins. Callables that
/// match neither trigger a `static_assert` at the call site.
///
/// Internally: the rank-2 gradient is materialized eagerly via
/// `diff::gradient(u)` at order 2 only when `f` consumes it; `f` is then
/// evaluated pointwise into a `Field<double>` integrand and reduced by
/// `func::integrate(integrand, tag)`. Higher-accuracy gradients
/// (`<Order = 4>`) and higher-rank derivatives (Laplacian /
/// biharmonic of vector fields) are not exposed in Phase 15; callers
/// who want either can pre-compute them and feed the materialized
/// auxiliary fields into a callable that takes them as captures.
///
/// \tparam F     Callable type; deduced from the argument.
/// \tparam Tag   Reduction-strategy tag --- must be `Pairwise` or `Kahan`.
/// \param  u     Input vector field.
/// \param  f     Functional integrand.
/// \param  tag   Reduction-strategy tag forwarded to `func::integrate`
///               (defaulted to `Pairwise{}`).
/// \returns The discrete functional value as a `double`.
/// \since 0.15.0
template <class F, class Tag = Pairwise>
inline double evaluate(const core::Field<core::Vec3d>& u, F&& f, Tag tag = {}) {
  static_assert(std::is_same_v<Tag, Pairwise> || std::is_same_v<Tag, Kahan>,
                "func::evaluate: Tag must be func::Pairwise or func::Kahan");

  const auto& grid = u.getGrid();
  const int Nx = grid.getNx();
  const int Ny = grid.getNy();
  const int Nz = grid.getNz();

  core::Field<double> integrand(grid);

  if constexpr (std::is_invocable_r_v<double, F&, const core::Vec3d&, const core::Mat3d&>) {
    const auto grad_u = diff::gradient(u);
    for (int k = 0; k < Nz; ++k) {
      for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
          integrand(i, j, k) = static_cast<double>(f(u(i, j, k), grad_u(i, j, k)));
        }
      }
    }
  } else if constexpr (std::is_invocable_r_v<double, F&, const core::Vec3d&>) {
    for (int k = 0; k < Nz; ++k) {
      for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
          integrand(i, j, k) = static_cast<double>(f(u(i, j, k)));
        }
      }
    }
  } else {
    static_assert(detail::deferred_false<F>::value,
                  "func::evaluate: f must be invocable as f(const Vec3d&) "
                  "or f(const Vec3d&, const Mat3d&), with a return type "
                  "convertible to double");
  }

  return integrate(integrand, tag);
}

}  // namespace gridcalc::func
