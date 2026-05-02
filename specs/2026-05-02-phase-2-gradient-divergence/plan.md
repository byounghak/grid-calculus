# Plan: Phase 2 — Gradient and divergence (scalar)

## Group 1 — Lift `Vec3d` into a shared aliases header

- Create `include/gridcalc/core/EigenAliases.hpp`. Move the `using Vec3d = Eigen::Vector3d;` declaration there with the `\since 0.1.0` brief preserved (Phase 1 already exposes the name, this is just a relocation).
- Update `include/gridcalc/core/Grid.hpp` to include the new header instead of declaring the alias inline. Public name `gridcalc::core::Vec3d` is unchanged, so Phase 1 callers (tests, downstream code) keep compiling.
- No new public symbols beyond the relocation.

**Commit message:** `refactor(core): lift Vec3d alias to core/EigenAliases.hpp`

## Group 2 — Add 1st-derivative central-difference stencil table

- Add `include/gridcalc/stencil/FirstDerivative.hpp` exposing `stencil::FirstDerivative<Order>`. Same primary-undefined + specialization pattern as `Coefficients<Order>`.
- Specialize `FirstDerivative<2>` with `radius = 1`, `offsets = {-1, 0, 1}`, `weights = {-0.5, 0.0, 0.5}` (consumer scales by `1/h`).
- Doxygen brief notes that `FirstDerivative<4>` is deferred to Phase 7 alongside `Coefficients<4>`.

**Commit message:** `feat(stencil): add FirstDerivative<Order=2> central-difference table`

## Group 3 — Implement `gradient` returning `Field<Vec3d>`

- Add `include/gridcalc/diff/Gradient.hpp` with `inline Field<Vec3d> gradient(const Field<double>&)`.
- Apply `FirstDerivative<2>` along x, y, z independently; pack per-point result as `Vec3d(dx/hx, dy/hy, dz/hz)`. Boundary wrap is delegated to the input `Field`'s `Policy` (`Periodic` at Phase 2).
- Doxygen block: `\brief`, `\param`, `\returns`, `\since 0.2.0`.

**Commit message:** `feat(diff): add 2nd-order periodic gradient on Field<double>`

## Group 4 — Implement `divergence` consuming `Field<Vec3d>`

- Add `include/gridcalc/diff/Divergence.hpp` with `inline Field<double> divergence(const Field<Vec3d>&)`.
- Same `FirstDerivative<2>` stencil; compute `∂Vx/∂x + ∂Vy/∂y + ∂Vz/∂z` pointwise.
- Doxygen block as above.

**Commit message:** `feat(diff): add 2nd-order periodic divergence on Field<Vec3d>`

## Group 5 — Convergence + round-trip identity tests

- Add `test/gradient_test.cpp`:
  - `GradientTest.RecoversTrigGradient` — at `N=32`, `ψ = sin(x) sin(y) sin(z)` on `L=2π`, expected ∇ψ has known analytical form. Relative max-norm error tolerance per the analytical leading-order bound.
  - `GradientTest.SecondOrderConvergence` — sweep `N ∈ {16, 32, 64}`, log-log slope of relative max-norm error in `[1.8, 2.2]`; halving-`h` ratio in `[3.5, 4.5]`.
  - `GradientTest.OutputGridMatchesInput` — shape preservation.
- Add `test/divergence_test.cpp`:
  - `DivergenceTest.RecoversAnalyticalDivergence` — at `N=32` on `V = (sin x cos y, cos x sin y, sin z)` (analytical `div V = 2 cos x cos y + cos z`). Pointwise tolerance.
  - `DivergenceTest.SecondOrderConvergence` — same sweep on the same V.
  - `DivergenceTest.RoundTripWithGradientMatchesLaplacian` — `divergence(gradient(ψ))` converges to the analytical `-(kx²+ky²+kz²) ψ` at order 2 over the same sweep. Note: `D₁∘D₁` and the direct 2nd-derivative stencil differ at finite `h` (different leading constants, both `O(h²)`), so the test compares to the *analytical* Laplacian, not the discrete one.
  - `DivergenceTest.OutputGridMatchesInput`.
- Wire the two new test sources into `test/CMakeLists.txt`.

**Commit message:** `test(diff): add gradient/divergence convergence + round-trip tests`

## Group 6 — Bookkeeping

- Update `STATUS.md`: bump "Last updated", move Phase 2 row to Done, add any decisions worth surfacing (e.g., `FirstDerivative<Order>` as a sibling stencil table to `Coefficients<Order>`; `core/EigenAliases.hpp` as the new home for vector typedefs).
- Update `CHANGELOG.md` (create if missing) under `0.2.0`.

**Commit message:** `docs: refresh STATUS for Phase 2 completion + 0.2.0 changelog entry`
