# Validation: Phase 2 — Gradient and divergence (scalar)

## Build

- [ ] CI green on Ubuntu GCC and Ubuntu Clang (Phase 21 will widen to Apple Clang + MSVC).
- [ ] Warnings clean under `-Wall -Wextra -Wpedantic -Wconversion` (Eigen still propagated as `SYSTEM`).
- [ ] `clang-format` and `clang-tidy` pass.
- [ ] `gridcalc` target remains INTERFACE (header-only) — no source files added.
- [ ] `core::Vec3d` continues to resolve from existing Phase 1 includes (Grid.hpp transitively pulls in EigenAliases.hpp).

## Tests

- [ ] `GradientTest.RecoversTrigGradient` passes — relative max-norm error within the documented analytical bound at `N = 32`.
- [ ] `GradientTest.SecondOrderConvergence` passes — log-log slope in `[1.8, 2.2]` and ratio between consecutive refinements in `[3.5, 4.5]` on `N ∈ {16, 32, 64}`.
- [ ] `GradientTest.OutputGridMatchesInput` passes.
- [ ] `DivergenceTest.RecoversAnalyticalDivergence` passes on `V = (sin x cos y, cos x sin y, sin z)` at `N = 32`.
- [ ] `DivergenceTest.SecondOrderConvergence` passes (same tolerance pattern as gradient) on the vector field above.
- [ ] `DivergenceTest.RoundTripWithGradientMatchesLaplacian` passes — `divergence(gradient(ψ))` converges at order 2 to the analytical `-(kx² + ky² + kz²) ψ` over `N ∈ {16, 32, 64}`.
- [ ] `DivergenceTest.OutputGridMatchesInput` passes.
- [ ] All Phase 1 tests still pass (`Vec3d` relocation must not break anything).

## Documentation

- [ ] Every new public symbol carries a Doxygen block including `\brief`, `\param`, `\returns`, `\since 0.2.0` (`gradient`, `divergence`, `FirstDerivative<Order>`, `FirstDerivative<2>`, `EigenAliases.hpp`'s `Vec3d` re-declaration retains its `\since 0.1.0`).
- [ ] `CHANGELOG.md` entry under `0.2.0` lists: gradient operator, divergence operator, `FirstDerivative<2>` stencil table, `core/EigenAliases.hpp` relocation.
- [ ] `STATUS.md` updated: "Last updated" bumped, Phase 2 row marked Done, decisions worth surfacing added (separate `FirstDerivative` table, EigenAliases.hpp).

## Performance

- N/A at Phase 2. Performance work is centralized in Phase 20.
