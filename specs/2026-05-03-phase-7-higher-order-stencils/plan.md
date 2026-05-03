# Plan: Phase 7 — Higher-order accuracy stencils (4th-order)

## Group 1 — Add `Coefficients<4>` and `FirstDerivative<4>` specializations

- Update `include/gridcalc/stencil/CentralDifference.hpp` to specialize `Coefficients<4>`:
  - `radius = 2`
  - `offsets = {-2, -1, 0, 1, 2}`
  - `weights = {-1.0/12, 4.0/3, -5.0/2, 4.0/3, -1.0/12}` (unscaled — consumer multiplies by `1/h²`)
  - `\since 0.7.0` Doxygen tag
- Update `include/gridcalc/stencil/FirstDerivative.hpp` to specialize `FirstDerivative<4>`:
  - `radius = 2`
  - `offsets = {-2, -1, 0, 1, 2}`
  - `weights = {1.0/12, -2.0/3, 0.0, 2.0/3, -1.0/12}` (unscaled — consumer multiplies by `1/h`)
  - `\since 0.7.0`
- The primary templates remain undefined; instantiating `Coefficients<3>` or any other unsupported `Order` continues to produce a compile error at the call site.

**Commit message:** `feat(stencil): add Coefficients<4> and FirstDerivative<4> specializations`

## Group 2 — Template `diff::laplacian` on `Order`

- Refactor `include/gridcalc/diff/Laplacian.hpp` to:
  - `template <int Order = 2> inline core::Field<double> laplacian(const core::Field<double>& field) noexcept`.
  - Inside, replace `using Coeffs = stencil::Coefficients<2>` with `using Coeffs = stencil::Coefficients<Order>`. The rest of the loop (per-axis stencil application + per-axis `1/h²` scaling) stays identical because it iterates `m` over `Coeffs::offsets`.
  - `\since 0.1.0` on the function (already had it); `\since 0.7.0` on the new `Order` parameter, called out in the Doxygen `\tparam` block.
- All Phase-1 callers writing `laplacian(field)` continue to compile via the default `Order = 2`. C++17 allows calling a function template whose only template argument has a default with the no-template-arg syntax.

**Commit message:** `refactor(diff): laplacian templated on Order (default 2)`

## Group 3 — Template `diff::gradient` and `diff::divergence` on `Order`

- Refactor `include/gridcalc/diff/Gradient.hpp` analogously:
  - `template <int Order = 2> inline core::Field<core::Vec3d> gradient(...)`.
  - `using Coeffs = stencil::FirstDerivative<Order>;`.
  - Loop body unchanged.
- Refactor `include/gridcalc/diff/Divergence.hpp` analogously:
  - `template <int Order = 2> inline core::Field<double> divergence(...)`.
  - `using Coeffs = stencil::FirstDerivative<Order>;`.
- All Phase-2 callers continue to compile via the `Order = 2` default.

**Commit message:** `refactor(diff): gradient and divergence templated on Order (default 2)`

## Group 4 — Tests for the new order

- Add `test/higher_order_test.cpp` containing:
  - `HigherOrderTest.LaplacianOrder4Sweep` — sweep `N ∈ {16, 32, 64}` on `ψ = sin(x) sin(y) sin(z)`. Log-log slope of relative max-norm error vs `h` in `[3.5, 4.5]`. Mirrors Phase 1's order-2 sweep but at order 4.
  - `HigherOrderTest.GradientOrder4Sweep` — same sweep, on `gradient<4>(ψ)`. Slope in `[3.5, 4.5]`.
  - `HigherOrderTest.DivergenceOrder4Sweep` — same sweep, on `divergence<4>(V)` with `V = (sin x cos y, cos x sin y, sin z)` (the Phase-2 vector field). Slope in `[3.5, 4.5]`.
  - `HigherOrderTest.LaplacianOrder4MoreAccurateThanOrder2` — at fixed `N = 32`, `laplacian<4>(ψ)` has ≥ 10× smaller relative max-norm error than `laplacian<2>(ψ)` against the analytic reference.
  - `HigherOrderTest.GradientOrder4MoreAccurateThanOrder2` — same comparison, on the gradient.
  - `HigherOrderTest.DivergenceOrder4MoreAccurateThanOrder2` — same comparison, on the divergence.
  - `HigherOrderTest.Coefficients4WeightsAreCorrect` — `static_assert` (or `EXPECT_DOUBLE_EQ`) that `Coefficients<4>::weights` matches the expected hard-coded values; same for `FirstDerivative<4>::weights`. Catches typos in the table at test time as well as compile time.
- Existing Phase 1–2 tests (`laplacian_test.cpp`, `gradient_test.cpp`, `divergence_test.cpp`) are not modified — they continue to exercise the order-2 default path.
- Wire `higher_order_test.cpp` into `test/CMakeLists.txt`.

**Commit message:** `test(diff): add order-4 convergence + order-2 comparison tests`

## Group 5 — Phase 7 doc-notes

- `docs/user-guide/notes/phase-7-higher-order-stencils.md` — user-facing summary: how to switch order on each operator (`laplacian<4>(field)`, `gradient<4>(field)`, `divergence<4>(field)`), when to pick order 4 vs the default order 2, the larger stencil radius (2 vs 1) means higher cost per grid point. A worked example showing the 10–100× accuracy improvement at the same `N`.
- `docs/developer-note/notes/phase-7-higher-order-stencils.md` — five-section structure per `specs/CLAUDE.md` step 4d:
  1. **Theory**: central-difference accuracy ladder; symmetry constraints (odd derivatives use anti-symmetric stencils, even use symmetric); how Cartesian-orthogonal grid lets the per-axis stencil structure stay independent of dimension and per-axis cell size.
  2. **Math derivation**: explicit Taylor-matching derivation of `Coefficients<4>::weights = {-1/12, 4/3, -5/2, 4/3, -1/12}` for `∂²/∂x²` and `FirstDerivative<4>::weights = {1/12, -2/3, 0, 2/3, -1/12}` for `∂/∂x`; truncation error `O(h^4)` instead of `O(h^2)`; how the per-axis-independent weights work on anisotropic per-axis spacing without modification.
  3. **Algorithm**: how the existing operator loops generalize to `Order` template parameter; per-grid-point cost is `2*radius+1` reads per axis (so 5 vs 3 for order 4 vs 2); periodic wrap is delegated to `Field::Policy` as before.
  4. **Design decisions**: compile-time `<Order>` over runtime parameter (cite `requirements.md`); hard-coded weights over Fornberg generation; no order-3 specialization (Fornberg would compute it; the project's symmetric central-difference convention is order 2, 4, 6, ... only).
  5. **References**: Fornberg 1988 (DOI 10.1090/...) on the general weight-generation algorithm — used as the *theoretical* derivation tool even though we hard-code; LeVeque (Chapter 1, table 1.1 of 1st-derivative coefficients); a permanent-URL reference for the standard 4th-order Laplacian stencil (Wikipedia or NIST DLMF).

**Commit message:** `docs(notes): add Phase 7 User Guide and Developer Note notes`

## Group 6 — Bookkeeping

- Bump `CMakeLists.txt` `project(gridcalc VERSION 0.6.0)` → `0.7.0`.
- Update `test/version_test.cpp` (`ReturnsZeroSixZero` → `ReturnsZeroSevenZero`).
- `CHANGELOG.md`: new `## 0.7.0 — Phase 7: 4th-order accuracy stencils` entry.
- `specs/STATUS.md`: bump "Last updated"; move Phase 7 row to Done; update "Where the project stands" + "Next action" (Phase 8 — higher-order spatial derivatives, 3rd & 4th).

**Commit message:** `chore: bump version to 0.7.0 and refresh STATUS for Phase 7`
