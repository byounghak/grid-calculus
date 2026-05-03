# Requirements: Phase 7 — Higher-order accuracy stencils (4th-order)

## Roadmap reference

> ## Phase 7 — Higher-order accuracy stencils (4th-order)
>
> - **Goal.** The accuracy-switching mechanism, validated on existing operators.
> - **Deliverables.**
>   - 4th-order central-difference coefficients for 1st and 2nd derivatives.
>   - `Order` template / runtime parameter on `gradient`, `laplacian`.
>   - Convergence tests for both orders.
> - **Acceptance.** Log-log convergence slopes match advertised orders within tolerance.

## Decisions made for this feature

- **`Order` is a compile-time template parameter on the operators.** `diff::laplacian<Order>(field)`, `diff::gradient<Order>(field)`, `diff::divergence<Order>(field)`. The default is `Order = 2`, so all Phase 1–2 callers (which write `laplacian(field)` etc.) keep compiling unchanged — C++17 resolves a function-template call with all template parameters defaulted to the default values. Compile-time dispatch matches the project convention (`stencil::Coefficients<Order>` is already templated this way) and produces a compile error for unsupported orders at the call site rather than silently producing wrong stencil weights or failing at runtime.
- **Hard-coded 4th-order weights, derivation in the Developer Note.** `Coefficients<4>::weights = {-1/12, 4/3, -5/2, 4/3, -1/12}` (for `∂²/∂x²`) and `FirstDerivative<4>::weights = {1/12, -2/3, 0, 2/3, -1/12}` (for `∂/∂x`), both unscaled. The consumer scales by `1/h²` or `1/h` per axis. Two specializations, ~10 LOC each, fitting the existing template pattern. No code generation. Constexpr Fornberg generation is deferred until Phase 8/11 actually need many more orders — at that point a compile-time generator becomes worth the implementation cost; right now hard-coding is leaner.
- **Anisotropic per-axis spacing works without modification.** The 4th-order weights are unitless and depend only on the integer index offsets `[−2, −1, 0, +1, +2]`. The 3D operator applies the same weight table along each Cartesian axis and scales by that axis's `1/h_a²` (Laplacian) or `1/h_a` (gradient/divergence) — exactly the same per-axis structure as Phase 1's order-2 Laplacian. Non-orthogonal Bravais grids (which would need cross-derivative terms) are out of scope for Phase 7 because Phase 1's `Grid` is Cartesian-orthogonal by contract; non-orthogonal lattices arrive in Phase 15.
- **All three currently-shipped differential operators get `Order`.** `laplacian` (Phase 1), `gradient` (Phase 2), `divergence` (Phase 2). Closes the order-switching API for everything that exists today; no operator-by-operator follow-up phases needed.
- **Test coverage: 3 convergence sweeps + 3 order-2-vs-order-4 comparisons + 1 weight-table verification.** The convergence sweep on each operator at order 4 (slope in `[3.5, 4.5]`, mirroring Phase 1–2's order-2 sweep tolerances) is the **roadmap acceptance**. The order-2 vs order-4 comparison at fixed `N = 32` (order 4 ≥ 10× more accurate) demonstrates the practical benefit. The weight-table verification (`EXPECT_DOUBLE_EQ` on `Coefficients<4>::weights` and `FirstDerivative<4>::weights`) catches typos in the hard-coded values.
- **Existing tests are not modified.** Phase 1's `LaplacianTest.*` and Phase 2's `GradientTest.*` / `DivergenceTest.*` continue to exercise the order-2 default path, which is the same code as before (just resolved through a template default). No regression risk.
- **No order-3 specialization.** Symmetric central-difference stencils have even truncation order (2, 4, 6, …). Order-3 would be an asymmetric upwind / off-center scheme — out of scope for this phase.
- **Version bump: `0.7.0`** on merge. New public symbol class (`Coefficients<4>`, `FirstDerivative<4>`); existing operators gain a backward-compatible template parameter. SemVer minor bump.

## Non-goals

- **No order ≥ 6 stencils.** Phase 8 considers higher-order spatial derivatives (3rd, 4th); Phase 11 considers up-to-4th-order *gradient* (not 6th-order accuracy). Order 6+ accuracy would need a Fornberg generator, deferred until any phase actually requires it.
- **No constexpr Fornberg algorithm.** Hard-coded weights for the two new specializations are simpler. Fornberg is the right tool when a phase needs many orders.
- **No order-3 specialization.** Symmetric central-difference stencils have even truncation order only.
- **No runtime-order parameter.** Compile-time only. A user-facing config wanting runtime order selection can write a small dispatch wrapper.
- **No cross-derivative stencils** (e.g., `∂²/∂x∂y`). Phase 6 territory (biharmonic + mixed partials) — wait, that's Phase 8 (higher-order spatial derivatives) per the roadmap. Either way not Phase 7.
- **No non-orthogonal-grid extension.** Phase 15.

## Deferred questions

- **Constexpr Fornberg generator.** Worth implementing when a phase introduces more than two new orders. Phase 8 (3rd, 4th derivatives) will be the natural decision point.
- **Order ≥ 6 specializations.** Out of scope until a phase needs them.
- **Mixed-partial / cross-derivative stencils.** Phase 8.
- **Runtime-order dispatch wrapper.** Add only if a downstream caller needs config-file-driven order selection.
