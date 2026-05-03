# Validation: Phase 7 — Higher-order accuracy stencils (4th-order)

## Build

- [ ] CI green on Ubuntu GCC and Ubuntu Clang.
- [ ] Warnings clean under `-Wall -Wextra -Wpedantic -Wconversion`.
- [ ] `clang-format` and `clang-tidy` pass (no new regressions vs the codebase-wide hygiene gap noted in Phase 2's STATUS).
- [ ] `gridcalc` target remains INTERFACE — no new source files.
- [ ] `Coefficients<3>` and any other unsupported `Order` continue to produce a compile error if instantiated (the primary template stays undefined; verified by inspection — no test exists because a successful compile of an unsupported instantiation would itself be the failure).

## Tests

- [ ] `HigherOrderTest.LaplacianOrder4Sweep` — log-log slope on `N ∈ {16, 32, 64}` lies in `[3.5, 4.5]`.
- [ ] `HigherOrderTest.GradientOrder4Sweep` — same.
- [ ] `HigherOrderTest.DivergenceOrder4Sweep` — same.
- [ ] `HigherOrderTest.LaplacianOrder4MoreAccurateThanOrder2` — at `N = 32`, `laplacian<4>` rel error is ≥ 10× smaller than `laplacian<2>` rel error.
- [ ] `HigherOrderTest.GradientOrder4MoreAccurateThanOrder2` — same.
- [ ] `HigherOrderTest.DivergenceOrder4MoreAccurateThanOrder2` — same.
- [ ] `HigherOrderTest.Coefficients4WeightsAreCorrect` — `Coefficients<4>::weights` and `FirstDerivative<4>::weights` match the expected hard-coded values (catches typos).
- [ ] All Phase 1–6 tests still pass without modification (they call `laplacian(field)` etc., which now resolve through the `Order = 2` default).

## Documentation

- [ ] Every new public symbol carries a Doxygen block. The `Coefficients<4>` and `FirstDerivative<4>` specializations get `\since 0.7.0`. The operators (`laplacian`, `gradient`, `divergence`) gain a `\tparam Order` line and a `\since 0.7.0` note alongside the existing function `\since` (`0.1.0` for `laplacian`, `0.2.0` for `gradient`/`divergence`).
- [ ] User Guide note `docs/user-guide/notes/phase-7-higher-order-stencils.md` covering: how to switch order on each operator, when to pick order 4 vs 2, larger stencil radius means higher per-point cost, a worked example demonstrating the 10–100× accuracy improvement.
- [ ] Developer Note `docs/developer-note/notes/phase-7-higher-order-stencils.md` with the five-section structure (Theory · Math derivation · Algorithm · Design decisions · References) per `specs/CLAUDE.md` step 4d. References section non-empty: at least one external citation.
- [ ] `CHANGELOG.md` entry under `0.7.0` lists the new stencil specializations and the `Order` template parameter on each operator.
- [ ] `STATUS.md` updated: Phase 7 row marked Done; Phase 8 stated as the new Next action.
- [ ] CI `render-docs` job produces fresh aggregate PDFs that include the Phase 7 chapter in both books.

## Performance

- N/A at Phase 7. The 4th-order stencil reads `2 * radius + 1 = 5` samples per axis instead of `3`, so per-point cost is ~`5/3 ≈ 1.67×` higher; documented in the Developer Note. Phase 20 owns optimization.
