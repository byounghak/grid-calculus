# Validation: Phase 5 — Explicit Euler diffusion solver

## Build

- [ ] CI green on Ubuntu GCC and Ubuntu Clang.
- [ ] Warnings clean under `-Wall -Wextra -Wpedantic -Wconversion`.
- [ ] `clang-format` and `clang-tidy` pass (no new regressions vs the codebase-wide hygiene gap noted in Phase 2's STATUS).
- [ ] `gridcalc` target remains INTERFACE — no new source files.
- [ ] No new top-level dependency added (no OpenMP, no MPI; explicit Euler is single-threaded).

## Tests

- [ ] `DiffusionTest.TrigEigenfunctionDecaysToAnalytic` passes at `N = 32, D = 0.1, dt = 0.05, n_steps = 100, T = 5.0`. Relative max-norm error vs `exp(-3 D T) ψ_0` is `< 5e-2`.
- [ ] `DiffusionTest.TrigEigenfunctionSecondOrderConvergence` passes — log-log slope on `N ∈ {16, 32, 64}` (with `dt ∝ h²`) lies in `[1.8, 2.5]`.
- [ ] `DiffusionTest.GaussianSanity` passes — Gaussian IC after 100 steps has finite values, conserves mass within `1e-10` relative (via `func::integrate`), and the peak max-norm has decreased.
- [ ] `DiffusionTest.CFLViolationThrows` passes — `diffuse` with violating parameters throws `std::invalid_argument`.
- [ ] `DiffusionTest.GenericExplicitEulerOnZeroRhs` passes — zero-RHS integration leaves `psi` bit-identical.
- [ ] `DiffusionTest.RejectsNegativeNSteps` passes — both `explicitEuler` and `diffuse` throw on `n_steps < 0`.
- [ ] All Phase 1–4 tests still pass.

## Documentation

- [ ] Every new public symbol (`solve::explicitEuler`, `solve::diffuse`) carries a Doxygen block with `\brief`, `\param`, `\returns`, `\throws`, `\since 0.5.0`. Internal helper `detail::checkExplicitDiffusionCFL` carries an internal-style comment but does not need full Doxygen.
- [ ] User Guide note at `docs/user-guide/notes/phase-5-explicit-euler.md` covering: signatures, CFL formula and throw behavior, a worked example.
- [ ] Developer Note at `docs/developer-note/notes/phase-5-explicit-euler.md` with the five-section structure (Theory · Math derivation · Algorithm · Design decisions · References) per `specs/CLAUDE.md` step 4d. References section non-empty: at least one external citation.
- [ ] `CHANGELOG.md` entry under `0.5.0` lists `solve::explicitEuler`, `solve::diffuse`, and the CFL-throw behavior.
- [ ] `STATUS.md` updated: Phase 5 row marked Done; Phase 6 stated as the new Next action.
- [ ] CI `render-docs` job produces fresh aggregate PDFs that include the Phase 5 chapters in both books.

## Performance

- N/A at Phase 5. Per-step allocation is documented as a known cost; Phase 20 owns the optimization.
