# Validation: Phase 11 — Functional with up-to-4th-order gradients

## Build

- [ ] CI green on Ubuntu GCC and Ubuntu Clang.
- [ ] `docs-build` CI job green (full LaTeX User Guide + Developer Note + Doxygen API reference).
- [ ] `docs-lint` CI job green (Doxygen `WARN_IF_UNDOCUMENTED` + `scripts/check-since.py`).
- [ ] Warnings clean under `-Wall -Wextra -Wpedantic -Wconversion`.
- [ ] `clang-format` and `clang-tidy` pass.
- [ ] `gridcalc` target remains INTERFACE — no new source files.
- [ ] No new top-level dependency added.
- [ ] Build green with `-DGRIDCALC_USE_FFT=OFF` (the spectral cross-check test must `GTEST_SKIP` cleanly).

## Tests

- [ ] `PfcFunctionalTest.HandComputedAtN64` passes — relative error vs $F_{\text{ref}} = 2.05546875\,\pi^3$ is `< 2e-2`. **(Roadmap acceptance bar.)**
- [ ] `PfcFunctionalTest.SecondOrderConvergence` passes — log-log slope in `[1.7, 2.3]`, halving ratio in `[3.5, 4.5]` over `N ∈ {16, 32, 64}`.
- [ ] `PfcFunctionalTest.SpectralCrossCheck` passes — $|F_{\text{spectral}} - F_{\text{ref}}| / |F_{\text{ref}}| < 1\text{e-}12$ at `N = 32`. (Skipped under `-DGRIDCALC_USE_FFT=OFF`.)
- [ ] `PfcFunctionalTest.PsiOnlyArityStillWorks` passes — `∫ ψ² dV ≈ π³` to spectral precision.
- [ ] `PfcFunctionalTest.GinzburgLandauStillWorks` passes — Phase 4's GL hand-computed value reproduced via the 3-arg dispatch path.
- [ ] All Phase 1–10 tests still pass without modification (no regression on the 209-test baseline).
- [ ] Compile-time error on a callable with no acceptable arity (verified manually; not part of the test suite). The static-assert message enumerates all four supported arities.

## Documentation

- [ ] `func::evaluate`'s function-level Doxygen block lists the new 4-arg signature `f(double, const Vec3d&, double, double)` and explicitly identifies the fourth argument as $\nabla^4\psi$ via `diff::biharmonic`.
- [ ] `func::evaluate`'s `\since` reads `0.4.0 (function); 0.11.0 (4-arg arity)`. The `\file` block reflects the same.
- [ ] User Guide chapter `docs/user-guide/chapters/11-higher-order-functional.tex` added, covering the new arity and the PFC worked example end-to-end.
- [ ] Developer Note chapter `docs/developer-note/chapters/10-higher-order-functional.tex` added with the five-section structure (Theory · Math derivation · Algorithm · Design decisions · References). References section non-empty: at least one external citation (Elder & Grant 2004, Phys. Rev. E 70, DOI link).
- [ ] `docs/user-guide/main.tex` `\input{}` order updated; `11-cahn-hilliard.tex` renamed to `12-cahn-hilliard.tex`; `12-diamond-lattice.tex` renamed to `13-diamond-lattice.tex`. Their `\label{chap:...}` anchors are preserved.
- [ ] `docs/developer-note/main.tex` `\input{}` list includes the new chapter 10.
- [ ] Title-page version strings on both LaTeX skeletons bumped from `0.10.0` to `0.11.0`.
- [ ] `CHANGELOG.md` entry under `0.11.0` lists the extended `func::evaluate` arity and the PFC test set.
- [ ] `STATUS.md` updated: Last-updated date, Phase 11 row marked Done, Phase 12 stated as the new Next action, the contracted-scalar-only and ET-deferral decisions captured under "Decisions worth knowing."

## Performance

- N/A at Phase 11. Performance work is centralized in Phase 20.
