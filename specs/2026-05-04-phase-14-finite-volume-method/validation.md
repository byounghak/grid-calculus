# Validation: Phase 14 — Finite volume method (FVM)

## Build

- [ ] CI green on Ubuntu GCC and Ubuntu Clang (Apple Clang + MSVC remain Phase 22 — renumbered).
- [ ] `docs-build` and `docs-lint` jobs green; new chapter 14 (User Guide) and chapter 13 (Developer Note) render with internal cross-references resolving.
- [ ] `\since`-tag lint passes — every new public symbol (`fvm::cellLaplacian`, the heterogeneous-D `solve::diffuse` overload) carries `\since 0.14.0`.
- [ ] Phase 9's FD–FFT cross-check fixture still passes — Phase 14 adds no new FD operators.
- [ ] Warnings clean under `-Wall -Wextra -Wpedantic -Wconversion`.
- [ ] `clang-format` / `clang-tidy` clean.

## Tests

### `fvm::cellLaplacian` (constant-D agreement)

- [ ] `AgreesWithFdLaplacianAtConstantUnitD` — `fvm::cellLaplacian(psi, D=1.0_field)` matches `diff::laplacian(psi)` to round-off on a periodic trig manufactured solution.
- [ ] `AgreesWithFdLaplacianAtConstantD` — same with `D = c_field` for `c = 0.7`; result equals `c · diff::laplacian(psi)`.
- [ ] `ConstantPsiGivesZeroRhs` — `psi = const` ⇒ RHS is zero (FVM is conservative).

### `fvm::cellLaplacian` (heterogeneous-D)

- [ ] `HarmonicMeanFaceCoefficient` — D-jump test: 2-cell face with `D_low` / `D_high` produces the harmonic-mean face flux to round-off.
- [ ] `AnalyticalSolutionConvergence_Order2` — `O(h²)` slope on `N ∈ {16, 32, 64}` for the manufactured solution `D(x) = 1 + 0.5 cos(x)`, `ψ(x) = sin(x)` (slope ∈ `[1.7, 2.3]`).
- [ ] `MassConservation_PeriodicSum` — `Σ fvm::cellLaplacian(ψ, D) = 0` to round-off (telescoping fluxes on a periodic domain).

### `solve::diffuse` heterogeneous-D overload

- [ ] `HeterogeneousD_AgreesWithConstantDOnUniformField` — running both the constant-D and the heterogeneous-D paths with the same uniform `D` produces the same trajectory after `n_steps` to round-off.
- [ ] `HeterogeneousD_RejectsCFLViolation` — `dt` above the heterogeneous CFL bound (`D_max · dt · Σ_a 1/h_a² > Tag::diffusionCFLLimit`) throws `std::invalid_argument`.
- [ ] `HeterogeneousD_MassConservedOverTimeStepping` — `Σ ψ` conserved to round-off across many integration steps.

### Acceptance — end-to-end manufactured-solution

- [ ] `HeterogeneousD_AnalyticalEndToEnd` — for a 1D-flavored manufactured solution with closed-form decay under heterogeneous diffusion, the discrete trajectory at `t = 1.0` matches the analytical answer to `O(h²)` tolerance on `N ∈ {32, 64}`. Discharges the roadmap acceptance bar.

### Regression

- [ ] All previously-passing tests remain green (243 prior on `clang-debug`, 167 prior on `clang-debug-nofft`). Phase 14 adds ~10 new tests across three files.

## Documentation

- [ ] User Guide chapter 14 (`docs/user-guide/chapters/14-finite-volume-method.tex`) — full chapter with motivation (FVM vs. FD), cell-flux discretization, harmonic-mean averaging, heterogeneous-`solve::diffuse` API, and the manufactured-solution acceptance test as a worked example.
- [ ] Developer Note chapter 13 (`docs/developer-note/chapters/13-finite-volume-method.tex`) — five-section structure with non-empty References (LeVeque 2007, Patankar 1980, Hundsdorfer & Verwer 2003, Eigen project docs at minimum).
- [ ] Both `main.tex` files updated with `\input{}` lines.
- [ ] Both LaTeX skeletons rebuild cleanly via `scripts/build-docs.sh`.
- [ ] `CHANGELOG.md` has a new `## 0.14.0 — Phase 14` block. The block explicitly calls out the roadmap renumber so future readers can correlate "Phase 14" in earlier blocks with the renamed phase.

## Roadmap renumber

- [ ] `specs/roadmap.md` — Phase 14 entry inserted; existing Phases 14–22 renumbered to 15–23; cross-references updated.
- [ ] `specs/STATUS.md` — "Next action" rewritten for the renumbered Phase 15; "Decisions worth knowing" cross-references updated; "Phase progress" table reflects the new numbering.
- [ ] `docs/user-guide/chapters/13-vector-tensor-fields.tex` — "Phase 14" → "Phase 15".
- [ ] `docs/developer-note/chapters/12-vector-tensor-fields.tex` — "Phase 14" → "Phase 15".
- [ ] `test/integrate_test.cpp` — comment "Phase 20" → "Phase 21".
- [ ] Sanity grep: no stale references to the pre-renumber phase numbers remain in the live narrative (CHANGELOG historical blocks intentionally untouched).

## Bookkeeping

- [ ] `project(... VERSION 0.14.0)` in root `CMakeLists.txt`.
- [ ] `test/version_test.cpp` asserts `0.14.0`.
- [ ] `specs/STATUS.md` updated: Phase 14 row → Done; "Last updated" date bumped; PR list extended.
- [ ] PR opened, every CI check green.
