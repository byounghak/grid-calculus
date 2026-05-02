# Validation: Phase 4 — Functional evaluation, callable API (scalar, periodic)

## Build

- [ ] CI green on Ubuntu GCC and Ubuntu Clang.
- [ ] Warnings clean under `-Wall -Wextra -Wpedantic -Wconversion`.
- [ ] `clang-format` and `clang-tidy` pass (no new regressions vs the codebase-wide hygiene gap noted in Phase 2's STATUS).
- [ ] `gridcalc` target remains INTERFACE — no new source files.
- [ ] No new top-level dependency added.

## Tests

- [ ] `FunctionalTest.GinzburgLandauHandComputed` passes at `N = 64` — relative error vs `(507/64)π³` is `< 1e-2`.
- [ ] `FunctionalTest.GinzburgLandauSecondOrderConvergence` passes — log-log slope in `[1.8, 2.2]`, ratio in `[3.5, 4.5]` over `N ∈ {16, 32, 64}`.
- [ ] `FunctionalTest.PsiOnlyArity` passes — `∫ ψ² dV ≈ π³` to spectral precision (1-arg dispatch path).
- [ ] `FunctionalTest.PsiAndGradArity` passes — `∫ (1/2)|∇ψ|² dV` matches `(3/2)π³` to `< 1e-2` relative (2-arg dispatch path).
- [ ] `FunctionalTest.PsiGradLapArity` passes — `∫ ψ·∇²ψ dV` matches `-3π³` to `< 1e-2` relative (3-arg dispatch path).
- [ ] `FunctionalTest.KahanTagAgrees` passes — Pairwise and Kahan results within `1e-13` relative on the GL setup.
- [ ] Compile-time error on a callable with no acceptable arity (verified manually; not part of the test suite).
- [ ] All Phase 1 + Phase 2 + Phase 3 tests still pass.

## Documentation

- [ ] Every new public symbol (`func::evaluate`) carries a Doxygen block with `\brief`, `\param`, `\returns`, `\since 0.4.0`.
- [ ] User Guide note added at `docs/user-guide/notes/phase-4-functional-evaluation.md`. Includes the GL worked example walked through end-to-end (the "tutorial doc page" deliverable).
- [ ] Developer Note added at `docs/developer-note/notes/phase-4-functional-evaluation.md` with the five-section structure (Theory · Math derivation · Algorithm · Design decisions · References) per `specs/CLAUDE.md` step 4d. References section non-empty: at least one external citation.
- [ ] `CHANGELOG.md` entry under `0.4.0` lists `func::evaluate`.
- [ ] `STATUS.md` updated: Phase 4 row marked Done; Phase 5 stated as the new Next action.

## Performance

- N/A at Phase 4. Performance work is centralized in Phase 20.
