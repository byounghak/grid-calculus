# Requirements: Phase 8 — Higher-order spatial derivatives (3rd, 4th)

## Roadmap reference

> ## Phase 8 — Higher-order spatial derivatives (3rd, 4th)
>
> - **Goal.** Mixed partials and the biharmonic operator.
> - **Deliverables.**
>   - `diff/MixedPartial.hpp` — $\partial^2/\partial x_i\partial x_j$.
>   - `diff/Biharmonic.hpp` — $\nabla^4$.
>   - 4th-order accuracy variants of both.
> - **Acceptance.** Convergence tests pass; biharmonic of $\sin(\mathbf{k}\cdot\mathbf{x})$ recovers $|\mathbf{k}|^4$.

## Decisions made for this feature

- **`specs/roadmap.md` Phase 8 entry edited in the same branch** to match the expanded scope below. New deliverables (`stencil/ThirdDerivative.hpp`, `stencil/FourthDerivative.hpp`, `diff/detail/MultiIndexDerivative.hpp`, `diff/ThirdOrder.hpp`, `diff/FourthOrder.hpp`) and the d-prefix naming convention are documented in the roadmap so future phases inherit the right contract; the original "mixed partials + biharmonic" wording was the bullet-point summary of a smaller scope decided before Phase 7's outcome was known.
- **Scope expanded to all unique scalar partials in 3D up to total rank 4.** The roadmap deliverables (rank-2 mixed partials + $\nabla^4$) are kept verbatim; on top of them this phase ships every distinct $\partial^{n_x+n_y+n_z}/\partial x^{n_x} \partial y^{n_y} \partial z^{n_z}$ with $1 \leq n_x+n_y+n_z \leq 4$ that is not already shipped (i.e., we do **not** re-ship rank 1 — that is `gradient` from Phase 2, axis-extracted via the existing `Vec3d`-valued result). Counts: 6 rank-2 (3 single-axis + 3 mixed), 10 rank-3, 15 rank-4. Total: **31 named partial-derivative functions** plus `biharmonic`. The expansion is the user's call after weighing the alternative of "single-axis only" — see Q1 of the Phase 8 spec round.
- **Naming convention: `d<rank>d<axis-with-exponent>...`, lexicographic axis order, exponent `1` omitted.** Examples: `d3dx2dy` is $\partial^3/\partial x^2 \partial y$; `d3dxdydz` is $\partial^3/\partial x \partial y \partial z$; `d4dx2dy2` is $\partial^4/\partial x^2 \partial y^2$; `d4dx4` is $\partial^4/\partial x^4$. Each function has signature
  ```cpp
  template <int Order = 2>
  core::Field<double> dXdY...(const core::Field<double>& psi);
  ```
  with `Order` defaulting to 2 (consistent with Phase 7's `laplacian<Order=2>` etc.) and an `Order = 4` specialization shipped in this phase. Unsupported orders trip a compile error inside the underlying `stencil::*Derivative<Order>` table — same mechanism as Phases 1, 2, 7.
- **Rank-2 partials use the same d-prefix naming** (`dxx, dyy, dzz, dxy, dxz, dyz`), not a runtime-axis `mixedPartial(psi, i, j)`. Internal consistency is more valuable than mirroring the roadmap's `diff/MixedPartial.hpp` filename literally — and the filename is preserved as the *home* for the six rank-2 functions.
- **Biharmonic = primary `biharmonic<Order>(psi)` + alias `d4<Order>(psi)`.** `biharmonic` is the physics-conventional name, kept as the primary entry point. `d4(psi)` is added as a thin alias that completes the d-prefix family at the contracted-rank level: $\nabla^4 = d^4 \psi$ (in the `d<rank>` shorthand with no axis suffix being interpreted as the contracted scalar operator on an even rank). Laplacian is **not** renamed and **no `d2` alias is added** — Phase 1's API stays put.
- **Implementation: direct fused tensor-product stencils.** Each multi-axis partial is computed in one pass over the grid by summing over the outer product of the relevant 1D weight tables. No intermediate `Field` allocations. The biharmonic uses a direct sum-of-six-terms stencil ($\partial_x^4 + \partial_y^4 + \partial_z^4 + 2\sum_{i<j}\partial_i^2\partial_j^2$) computed in a single pass, not `laplacian(laplacian(psi))`. The rejected alternative (composition through an intermediate `Field`) inherits the stride-2 leading-error behavior already documented in `STATUS.md` for `divergence(gradient(psi))` — same convergence order but different leading constant; we avoid that drift by going direct.
- **New stencil weight tables: `stencil::ThirdDerivative<Order>` and `stencil::FourthDerivative<Order>`.** Sibling templates to `stencil::Coefficients<Order>` (= 2nd derivative, Phase 1) and `stencil::FirstDerivative<Order>` (Phase 2). Each is hard-coded — no Fornberg generator. New file pair: `include/gridcalc/stencil/ThirdDerivative.hpp` and `include/gridcalc/stencil/FourthDerivative.hpp`. The unification meta-decision (a single `Derivative<Rank, Order>` 2-parameter template) was considered for the second time in this phase and again deferred — keeping the existing pattern is cheaper and the project still has no use case where a `Rank` template parameter would be passed dynamically.
- **Hard-coded weights, no Fornberg generator.** Same call as Phase 7. Eight new specializations (`ThirdDerivative<2>`, `ThirdDerivative<4>`, `FourthDerivative<2>`, `FourthDerivative<4>` × 1 each derivative table) totaling ~30 LOC. Fornberg becomes worth its weight when a phase needs orders ≥ 6 or many ranks at once; today's set is small enough that a typo audit by `EXPECT_DOUBLE_EQ` against the published Fornberg-derived values is the lighter approach.
- **Convergence-sweep grid sizes: `N ∈ {16, 32, 64}`.** Same triplet as Phases 1, 2, 7 — keeps the slope window familiar and CI duration manageable. The order-4 sweep on the most stencil-intensive operator (the rank-4 `d4dx2dy2<4>` reaches `2 * radius + 1 = 5` samples per axis, so a 5×5 = 25-point cross stencil) at `N = 16` retains a halo-to-grid ratio of `2/16 = 12.5 %` — safe under `IndexPolicy::Periodic`.
- **Test coverage: 1 convergence sweep per (operator, Order) pair + biharmonic acceptance + 8 weight-table audits.** That's 31 partials × 2 orders + biharmonic × 2 orders = **64 convergence sweeps**, each on the trigonometric manufactured solution $\psi = \sin(x)\sin(y)\sin(z)$ over the periodic $[0, 2\pi)^3$ box (analytical reference computable in closed form for every partial — see the Developer Note). Tests are bundled into a parameterized `HigherOrderDerivativeTest` GoogleTest fixture that loops over operator-tags rather than 64 individual `TEST` macros. Plus the explicit roadmap acceptance: `BiharmonicTest.RecoversKToTheFourth` on $\psi = \sin(\mathbf{k}\cdot\mathbf{x})$ for two non-trivial $\mathbf{k}$ vectors. Plus 8 weight-table verifications (`EXPECT_DOUBLE_EQ` on `ThirdDerivative<2>`, `ThirdDerivative<4>`, `FourthDerivative<2>`, `FourthDerivative<4>` weights).
- **Existing tests are not modified.** Phase 1–7 tests continue to exercise `laplacian`, `gradient`, `divergence` at their default `Order = 2` paths. Phase 8 adds new files only.
- **Version bump: `0.8.0`.** New public symbols (32 free functions in `gridcalc::diff` + two new stencil templates with two specializations each + `d4` alias). No breaking change to any prior-phase API. SemVer minor bump.

## Non-goals

- **No spectral verification.** That is Phase 9. Phase 8 ships a self-consistent FD test suite (analytical reference per partial); the FFT cross-check joins on the next phase.
- **No order ≥ 6 stencils.** Same call as Phase 7; needs a Fornberg generator the project has not yet built.
- **No order-3 specializations.** Symmetric central-difference stencils have even truncation order (2, 4, 6, …). Order 3 is asymmetric / upwind territory, out of scope.
- **No vector- or tensor-valued partial-derivative entry points.** Phase 11 (up-to-4th-order *gradient*) and Phase 13/14 (vector / tensor fields) own those. Phase 8's `d<rank>...` family is `Field<double> → Field<double>` only.
- **No constexpr Fornberg generator.** Deferred again. The next natural decision point is whichever phase first needs many new (rank, order) combinations — likely Phase 11.
- **No non-orthogonal-grid extension.** Phase 15.
- **No runtime-order or runtime-multi-index dispatch.** All naming and order selection is compile-time. A wrapper layer can be added downstream by anyone who needs config-file-driven dispatch.

## Deferred questions

- **Constexpr Fornberg generator.** Deferred to Phase 11 unless an earlier phase changes the calculus.
- **Order-3 / odd-order asymmetric stencils.** Out of scope until a phase has a use case (probably never on periodic-only grids).
- **Vector- and tensor-valued higher-order partials** (gradient of rank ≥ 2 input data). Phase 11 / 13 / 14.
- **Spectral cross-check fixture extending to Phase 8 operators.** Phase 9 — the parameterized FD–FFT discrepancy test in Phase 9 is required to cover every Phase 8 operator at every advertised accuracy order, per the roadmap's "permanent CI gate" wording.
