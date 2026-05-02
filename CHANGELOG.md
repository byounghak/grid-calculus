# Changelog

All notable changes to gridcalc are documented in this file. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project adheres to [Semantic Versioning](https://semver.org/).

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
