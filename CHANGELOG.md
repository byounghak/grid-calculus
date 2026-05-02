# Changelog

All notable changes to gridcalc are documented in this file. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project adheres to [Semantic Versioning](https://semver.org/).

## 0.5.0 — Phase 5: explicit Euler diffusion solver

### Added

- `gridcalc::solve::explicitEuler<Rhs>(Field<double>&, Rhs&&, double dt, int n_steps)`
  — generic forward-Euler integrator; advances
  $\psi^{n+1} = \psi^n + dt \cdot \mathrm{rhs}(\psi^n)$ for an arbitrary
  RHS callable. Mutates `psi` in place. SFINAE-validates the callable
  signature; throws `std::invalid_argument` on negative `dt` or
  `n_steps`.
- `gridcalc::solve::diffuse(Field<double>&, double D, double dt, int n_steps)`
  — convenience driver for $\partial_t\psi = D\nabla^2\psi$. Internally
  calls `explicitEuler` with `diff::laplacian` as the RHS and `D * dt`
  as the effective step. Performs a von Neumann CFL check before
  integration:
  $D \cdot dt \cdot \sum_a (1/h_a^2) \le 1/2$.
  Throws `std::invalid_argument` on violation, with a message naming
  the offending values.
- New public namespace `gridcalc::solve` and directory
  `include/gridcalc/solve/`.

### Tests

- Trig eigenfunction `ψ_0 = sin(x)sin(y)sin(z)` decays to the analytical
  `exp(-3 D T) ψ_0` after 100 explicit-Euler steps; relative max-norm
  error `< 5e-2` at `N = 32`.
- 2nd-order convergence sweep on the same trig IC over `N ∈ {16, 32, 64}`
  with `dt ∝ h²`; log-log slope in `[1.8, 2.5]`.
- Gaussian sanity test: finite values, mass conserved (via
  `func::integrate`), peak max-norm decreased after integration.
- CFL violation throws.
- Generic `explicitEuler` with zero RHS leaves `psi` bit-identical.
- Negative `n_steps` throws on both entry points.

### Documentation

- New User Guide note `docs/user-guide/notes/phase-5-explicit-euler.md`
  with worked examples for `diffuse` and a custom-RHS use of
  `explicitEuler`.
- New Developer Note `docs/developer-note/notes/phase-5-explicit-euler.md`
  with the von Neumann derivation of the CFL bound and four external
  references (Strikwerda DOI, LeVeque DOI, Wikipedia heat-equation
  permanent URL, Numerical Recipes ISBN).

### Changed

- `specs/CLAUDE.md` step 4a: new "Commit-message rules" subsection
  forbidding `Co-Authored-By:` trailers in this repository.

## 0.4.0 — Phase 4: functional evaluation, callable API (scalar, periodic)

### Added

- `gridcalc::func::evaluate(const Field<double>&, F&&, Tag = {}) -> double`
  — discrete functional $F[\psi] = \int f(\psi, \nabla\psi, \nabla^2\psi)\, dV$.
  The callable `f` is auto-detected at compile time as one of three
  arities:
  - `f(double psi)` — only `ψ`; the gradient and Laplacian are skipped.
  - `f(double psi, const Vec3d& grad)` — `ψ` and `∇ψ`; the Laplacian is
    skipped.
  - `f(double psi, const Vec3d& grad, double lap)` — full triplet.

  The reduction tag (`Pairwise` or `Kahan`) is forwarded to
  `func::integrate`. Implementation is eager: required derivatives are
  materialized via `diff::gradient` / `diff::laplacian` and a
  per-grid-point integrand `Field<double>` is reduced.

### Tests

- Ginzburg–Landau hand-computed reference at `N = 64` matches
  $F = (507/64)\pi^3 \approx 245.62$ within `1e-2` relative.
- 2nd-order convergence sweep on the GL functional over `N ∈ {16, 32, 64}`.
- Per-arity dispatch tests: `f(ψ) = ψ²`, `f(ψ, ∇ψ) = (1/2)|∇ψ|²`,
  `f(ψ, _, ∇²ψ) = ψ ∇²ψ`, all matching their analytical values.
- Pairwise and Kahan reducers agree on the GL setup within `1e-13` relative.

### Documentation

- New User Guide note `docs/user-guide/notes/phase-4-functional-evaluation.md`
  walking through the Ginzburg–Landau worked example end-to-end (the
  "tutorial doc page" deliverable from `roadmap.md` Phase 4).
- New Developer Note `docs/developer-note/notes/phase-4-functional-evaluation.md`
  with the five-section structure: theory of functional convergence,
  full closed-form derivation of the GL reference value, arity-dispatch
  algorithm, design decisions, and five external references.

## 0.3.0 — Phase 3: domain integration (∫ over the grid)

### Added

- `gridcalc::func::integrate(const Field<double>&, Tag = {}) -> double` —
  discrete periodic-midpoint integral
  $I_h = h_x h_y h_z \sum_{i,j,k} f_{ijk}$. Two reduction tags:
  - `gridcalc::func::Pairwise` (default): pairwise recursive summation,
    rounding error $O(\varepsilon \log n)$, 64-element base case.
  - `gridcalc::func::Kahan`: Neumaier-style compensated summation,
    rounding error independent of $n$ to first order.
- New public namespace `gridcalc::func` and corresponding directory
  `include/gridcalc/func/`.

### Tests

- Unit-field-recovers-domain-volume on an anisotropic non-cubic Grid.
- $\int \sin(kx)\, dV = 0$ over `[0, 2π]³`.
- Manufactured solution `cos(x) sin(2y) sin(3z)` integrates to 0.
- Pairwise and Kahan agree within `1e-13` relative on a non-trivial-scale
  integrand.
- Bit-equality across repeated invocations (placeholder for the
  `OMP_NUM_THREADS`-varying test that lands in Phase 20).

### Documentation

- New User Guide note `docs/user-guide/notes/phase-3-domain-integration.md`.
- New Developer Note `docs/developer-note/notes/phase-3-domain-integration.md`,
  first to use the five-section + external-citation rule introduced in
  `specs/CLAUDE.md` step 4d (commits `aabfbf6`, `2309ed4`).

## 0.2.0 — Phase 2: gradient and divergence (scalar)

### Added

- `gridcalc::diff::gradient(const Field<double>&) -> Field<Vec3d>` — 2nd-order
  central-difference gradient on a periodic scalar field. Returns Cartesian
  triple `(∂f/∂x, ∂f/∂y, ∂f/∂z)` per grid point.
- `gridcalc::diff::divergence(const Field<Vec3d>&) -> Field<double>` —
  2nd-order central-difference divergence on a periodic vector field.
- `gridcalc::stencil::FirstDerivative<Order>` — sibling to
  `Coefficients<Order>` (which is the second-derivative table). `Order = 2`
  specialization with weights `{-0.5, 0, 0.5}` (consumer scales by `1/h`).
  Primary template intentionally undefined; `Order = 4` deferred to Phase 7.
- `gridcalc/core/EigenAliases.hpp` — new home for project-local Eigen
  typedefs (`Vec3d` today; `Mat3d` and tensor aliases will land here when
  Phase 13/14 needs them).

### Changed

- `Vec3d` alias relocated from `core/Grid.hpp` to `core/EigenAliases.hpp`.
  Fully-qualified name `gridcalc::core::Vec3d` is unchanged; Phase 1 callers
  recompile without source changes.

### Tests

- Convergence-order tests pass for `gradient` and `divergence` on trig inputs
  over `N ∈ {16, 32, 64}` (log-log slope in `[1.8, 2.2]`, ratio in
  `[3.5, 4.5]`).
- Round-trip identity test: `divergence(gradient(ψ))` converges at order 2
  to the analytical Laplacian `-(kx²+ky²+kz²) ψ` on the same sweep.

## 0.1.0 — Phase 1: periodic 3D scalar grid + Laplacian

### Added

- `gridcalc::core::Grid` — Cartesian-orthogonal sampling grid with per-axis
  cell sizes; i-fastest storage layout.
- `gridcalc::core::Field<T, Policy = Periodic>` — 3D field on a `Grid` with a
  pluggable out-of-range index policy. Three constructors (zero-init,
  broadcast value, callable sampler).
- `gridcalc::core::IndexPolicy::Periodic` — wrap-around index policy.
- `gridcalc::stencil::Coefficients<Order>` (second-derivative table),
  `Order = 2` specialized; `Order = 4` deferred to Phase 7.
- `gridcalc::diff::laplacian` — 2nd-order central-difference periodic
  Laplacian.

### Tests

- Trig-eigenvalue recovery and 2nd-order convergence pass on `N ∈ {16, 32, 64}`.

## 0.0.1 — Phase 0: project scaffold

### Added

- CMake 3.20+ project skeleton, Eigen 3.4.0 + GoogleTest 1.14.0 pinned,
  `.clang-format` / `.clang-tidy` / `CMakePresets.json`, GitHub Actions CI on
  Ubuntu (GCC + Clang). Proprietary `LICENSE` (all rights reserved).
