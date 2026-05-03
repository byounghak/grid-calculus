# Changelog

All notable changes to gridcalc are documented in this file. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project adheres to [Semantic Versioning](https://semver.org/).

## 0.8.0 — Phase 8: Higher-order spatial derivatives (3rd, 4th)

### Added

- `gridcalc::stencil::ThirdDerivative<Order>` — central-difference weight
  table for $\partial^3/\partial x^3$. Order 2: $\{-1/2, 1, 0, -1, 1/2\}$
  (radius 2, anti-symmetric). Order 4: $\{1/8, -1, 13/8, 0, -13/8, 1, -1/8\}$
  (radius 3, anti-symmetric).
- `gridcalc::stencil::FourthDerivative<Order>` — central-difference weight
  table for $\partial^4/\partial x^4$. Order 2: $\{1, -4, 6, -4, 1\}$
  (radius 2, symmetric). Order 4:
  $\{-1/6, 2, -13/2, 28/3, -13/2, 2, -1/6\}$ (radius 3, symmetric).
- `gridcalc::diff::detail::mixedDerivative<Nx, Ny, Nz, Order>` — internal
  tensor-product helper that computes the multi-axis partial
  $\partial^{Nx+Ny+Nz}\psi / \partial x^{Nx}\partial y^{Ny}\partial z^{Nz}$
  in a single fused pass via the outer product of three 1D weight tables.
- **31 named partial-derivative functions** across three new headers,
  each templated on `<int Order = 2>` (orders 2 and 4 specialized):
  - `gridcalc/diff/MixedPartial.hpp` — rank-2 partials: `dxx`, `dyy`,
    `dzz`, `dxy`, `dxz`, `dyz` (6 functions). Roadmap deliverable
    `diff/MixedPartial.hpp` shipped under the d-prefix naming
    convention.
  - `gridcalc/diff/ThirdOrder.hpp` — rank-3 partials (10 functions):
    `d3dx3`, `d3dy3`, `d3dz3`, `d3dx2dy`, `d3dxdy2`, `d3dx2dz`,
    `d3dxdz2`, `d3dy2dz`, `d3dydz2`, `d3dxdydz`.
  - `gridcalc/diff/FourthOrder.hpp` — rank-4 partials (15 functions):
    `d4dx4`, `d4dy4`, `d4dz4`, `d4dx3dy`, `d4dxdy3`, `d4dx3dz`,
    `d4dxdz3`, `d4dy3dz`, `d4dydz3`, `d4dx2dy2`, `d4dx2dz2`,
    `d4dy2dz2`, `d4dx2dydz`, `d4dxdy2dz`, `d4dxdydz2`.
- `gridcalc::diff::biharmonic<Order>(psi)` — $\nabla^4\psi$ on a periodic
  scalar field, implemented as a fused six-term direct stencil
  ($\partial_x^4 + \partial_y^4 + \partial_z^4 +
  2\sum_{a<b}\partial_a^2\partial_b^2$) — not via
  `laplacian(laplacian(psi))`. Roadmap deliverable
  `diff/Biharmonic.hpp`.
- `gridcalc::diff::d4<Order>(psi)` — alias forwarding to `biharmonic`,
  completing the d-prefix family at the contracted-rank level.

### Tests

- **Weight-table audits**: `EXPECT_DOUBLE_EQ` on every weight in
  `ThirdDerivative<{2, 4}>` and `FourthDerivative<{2, 4}>` (catches typos
  in the hard-coded values).
- **62 convergence sweeps** across all 31 named partials × `Order ∈ {2, 4}`
  on the trigonometric manufactured solution $\psi = \sin(x)\sin(y)\sin(z)$
  over the periodic $[0, 2\pi)^3$ box: log-log slope on $N \in \{16, 32, 64\}$
  in $[\text{Order} - 0.5, \text{Order} + 0.5]$. Satisfies the roadmap's
  acceptance.
- **Biharmonic convergence sweeps** (orders 2 and 4) and **four explicit
  $\nabla^4\psi = |\mathbf{k}|^4\psi$ acceptance tests** at $N = 64$ on
  $\mathbf{k} = (1, 1, 1)$ and $(1, 2, 1)$ — satisfies the roadmap's
  "biharmonic of $\sin(\mathbf{k}\cdot\mathbf{x})$ recovers $|\mathbf{k}|^4$"
  bar.
- **Alias-parity tests** verifying `d4<Order>(psi) == biharmonic<Order>(psi)`
  element-wise.
- All Phase 1–7 tests continue to pass without modification.

### Documentation

- New User Guide note
  `docs/user-guide/notes/phase-8-higher-order-derivatives.md` covering the
  full naming convention, the 31 named functions plus `biharmonic` / `d4`,
  the calling syntax with `Order`, and a worked $\nabla^4$ example with the
  $|\mathbf{k}|^4$ spectral expectation cross-check.
- New Developer Note
  `docs/developer-note/notes/phase-8-higher-order-derivatives.md` with the
  five-section structure (Theory · Math derivation · Algorithm · Design
  decisions · References), explicit Taylor-matching derivation of the four
  new 1D weight tables, the biharmonic six-term decomposition, and three
  external references (Fornberg 1988 DOI, LeVeque 2007 DOI, NIST DLMF
  permanent URL).
- `specs/roadmap.md` Phase 8 entry rewritten to match the expanded scope
  decided in the spec round.

## 0.7.0 — Phase 7: 4th-order accuracy stencils

### Added

- `gridcalc::stencil::Coefficients<4>` — 4th-order central-difference
  weights for $\partial^2/\partial x^2$:
  $\{-1/12, 4/3, -5/2, 4/3, -1/12\}$ (radius 2). Truncation error
  $O(h^4)$ vs `Coefficients<2>`'s $O(h^2)$.
- `gridcalc::stencil::FirstDerivative<4>` — 4th-order central-difference
  weights for $\partial/\partial x$:
  $\{1/12, -2/3, 0, 2/3, -1/12\}$ (radius 2). Truncation error
  $O(h^4)$ vs `FirstDerivative<2>`'s $O(h^2)$.

### Changed

- `gridcalc::diff::laplacian`, `gradient`, and `divergence` are now
  templated on an `Order` parameter, defaulting to `2`. All Phase 1–6
  callers (which write `laplacian(field)`, etc.) keep compiling
  unchanged via the C++17 default-template-argument rule. Phase 7
  callers can write `laplacian<4>(field)`, `gradient<4>(field)`, or
  `divergence<4>(vfield)` to opt into the 4th-order stencil.

### Tests

- Order-4 convergence sweep on each operator over `N ∈ {16, 32, 64}`:
  log-log slope in `[3.5, 4.5]`. Satisfies the roadmap's "log-log
  convergence slopes match advertised orders" acceptance.
- Order-2-vs-order-4 absolute-accuracy comparison on each operator at
  `N = 32`: order 4 is at least 10× more accurate than order 2.
- Weight-table verification for `Coefficients<4>` and
  `FirstDerivative<4>` (catches typos in the hard-coded values).

### Documentation

- New User Guide note `docs/user-guide/notes/phase-7-higher-order-stencils.md`
  with the API switch, cost / accuracy table, and a worked example.
- New Developer Note `docs/developer-note/notes/phase-7-higher-order-stencils.md`
  with the explicit Taylor / moment-matching derivation of both 4th-order
  weight tables, the truncation-error analysis, and four external
  references (Fornberg 1988, LeVeque, NIST DLMF, Trefethen).

## 0.6.0 — Phase 6: RK4 + generic time integrator interface

### Added

- `gridcalc::solve::ExplicitEuler` and `gridcalc::solve::RK4` empty
  tag structs, each carrying a `static constexpr double diffusionCFLLimit`
  (Euler `0.5`, RK4 `0.6963`).
- `gridcalc::solve::integrate(psi, rhs, dt, n_steps, IntegratorTag{})` —
  generic tag-dispatched time integrator. Two overloads, one per tag.
- Classic 4-stage Runge–Kutta implementation in
  `solve/RK4.hpp`; per-step allocation is `~6` fields. Phase 20 will
  optimize via persistent scratch.
- `solve/Integrator.hpp` — one-stop aggregator that pulls in both
  `ExplicitEuler.hpp` and `RK4.hpp`.

### Changed

- `gridcalc::solve::diffuse` is now templated on the integrator tag:
  `diffuse(psi, D, dt, n_steps, Integrator{})`. The default tag is
  `ExplicitEuler{}`, so Phase 5 callers (no integrator argument) keep
  working unchanged. The CFL check now reads
  `Integrator::diffusionCFLLimit`, so RK4 callers get the larger
  stability bound automatically.
- `gridcalc::solve::explicitEuler` is now a thin wrapper that forwards
  to `solve::integrate(..., ExplicitEuler{})`. Phase 5 callers and
  tests are unaffected.
- The Phase-5-specific `solve::detail::checkExplicitDiffusionCFL` is
  superseded by the templated `solve::detail::checkDiffusionCFL<Tag>`.

### Tests

- RK4 matches the analytical heat-equation decay within `1e-2` relative
  at `N=32, T=5.0` (10× tighter than Phase 5's Euler at the same
  parameters).
- Combined refinement (`N ∈ {16, 32, 64}` with each integrator's own
  CFL-limited dt): both Euler and RK4 show log-log slope in `[1.8, 2.5]`
  (combined order is spatial-dominated at CFL-limited dt).
- Pure linear-decay ODE comparison: RK4's error is at least 10× smaller
  than Euler's at the same dt (isolates time-error from spatial-error).
- **Pure linear-decay ODE convergence sweep over `dt`**: RK4 log-log
  slope in `[3.5, 5.0]`, ratio `>= 8` per halving. Satisfies the
  roadmap's `error scales as Δt^4` acceptance.
- RK4 succeeds in the gap between Euler's CFL (`0.5`) and RK4's
  (`0.6963`) where Euler throws.
- Bit-identical output between Phase 5's `solve::explicitEuler(...)`
  wrapper and `solve::integrate(..., ExplicitEuler{})`.

### Documentation

- New User Guide note `docs/user-guide/notes/phase-6-rk4-integrator.md`.
- New Developer Note `docs/developer-note/notes/phase-6-rk4-integrator.md`
  with the four-stage Taylor-match derivation, the RK4 stability-region
  intersection at `~-2.7853`, and four external references
  (Hairer–Nørsett–Wanner, Butcher, Wikipedia RK methods, LeVeque).

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
