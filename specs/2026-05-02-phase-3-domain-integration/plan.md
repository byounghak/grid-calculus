# Plan: Phase 3 — Domain integration (∫ over the grid)

## Group 1 — Add `func::integrate` with pairwise default + Kahan overload

- Create `include/gridcalc/func/Integrate.hpp` (new `gridcalc::func` namespace; new `func/` directory under the existing public-include tree).
- Define empty tag structs `gridcalc::func::Pairwise` (default) and `gridcalc::func::Kahan` for reduction-strategy dispatch.
- Implement `inline double integrate(const Field<double>& field, Pairwise = {})` using a pairwise recursive sum (small-block base case ~64 elements, naive accumulator inside). Multiply the reduced sum by `grid.getCellVolume()` for the periodic midpoint/trapezoidal-on-equispaced-samples result.
- Implement `inline double integrate(const Field<double>& field, Kahan)` using Neumaier-style compensated summation (single pass, per-element compensation).
- No source files are added; `gridcalc` stays INTERFACE.
- Doxygen blocks on every public symbol (`integrate` overloads, `Pairwise`, `Kahan`) with `\since 0.3.0`.

**Commit message:** `feat(func): add integrate() with pairwise default + Kahan overload`

## Group 2 — Tests for integrate()

- Add `test/integrate_test.cpp`:
  - `IntegrateTest.UnitFieldRecoversDomainVolume` — `∫ 1 dV` equals `Lx*Ly*Lz` to machine precision on a non-cubic, anisotropic grid.
  - `IntegrateTest.SinOverPeriodIsZero` — `∫ sin(kx) dV` = 0 over `[0, 2π]³` for integer `kx`. Tolerance set by the periodic trapezoidal rule's spectral accuracy (≤ 1e-12 at moderate N).
  - `IntegrateTest.ManufacturedSolutionPeriodicProduct` — `∫ cos(x) sin(2y) sin(3z) dV` over `[0, 2π]³` equals 0 to machine precision; each axis integrates to zero independently. This is the manufactured-solution accuracy test.
  - `IntegrateTest.PairwiseAndKahanAgree` — for a typical integrand sized `64³`, `integrate(f, Pairwise{})` and `integrate(f, Kahan{})` agree within a tight relative tolerance (~1e-13).
  - `IntegrateTest.DeterministicAcrossInvocations` — calling `integrate()` twice on the same field returns bit-identical results. Placeholder for the cross-thread-count determinism contract; trivially satisfied at Phase 3 because the implementation is single-threaded. Phase 20 expands this to vary `OMP_NUM_THREADS`.
- Wire `integrate_test.cpp` into `test/CMakeLists.txt`.

**Commit message:** `test(func): add integrate accuracy + determinism tests`

## Group 3 — Phase 3 doc-notes (User Guide + Developer Note)

- `docs/user-guide/notes/phase-3-domain-integration.md` — user-facing API summary: signature of `integrate()`, semantics of `Pairwise` vs. `Kahan`, the periodic-midpoint quadrature rule, a worked example computing the volume of a 2π cube and the integral of a periodic trig product.
- `docs/developer-note/notes/phase-3-domain-integration.md` — full five-section structure per `specs/CLAUDE.md` step 4d:
  1. **Theory**: Riemann integration on a Cartesian periodic domain; equivalence of midpoint and trapezoidal rules on equispaced periodic samples; spectral accuracy of the periodic trapezoidal rule for analytic periodic integrands.
  2. **Math derivation**: derivation of the discrete formula `I_h = h_x h_y h_z ∑_{i,j,k} f_{ijk}`; exponential convergence for analytic periodic `f` (Trefethen–Weideman); error bounds for non-periodic `f` (Euler–Maclaurin) — flagged as out-of-scope until Phase 18 (Dirichlet/Neumann); rounding-error analysis showing pairwise has `O(ε log n)` rounding error vs naive `O(ε n)`.
  3. **Algorithm**: pairwise recursion structure, base-case threshold rationale; Neumaier compensation step.
  4. **Design decisions**: why both Pairwise and Kahan are exposed (cite `requirements.md`); why Phase 3 is single-threaded (cite Phase 20 follow-up).
  5. **References**: Trefethen & Weideman on the periodic trapezoidal rule (DOI); Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed. (textbook citation, §4.3 pairwise summation, §4.4 Kahan/Neumaier); Neumaier 1974 short note (DOI or full bibliographic entry).

**Commit message:** `docs(notes): add Phase 3 User Guide and Developer Note notes`

## Group 4 — Bookkeeping

- Bump `CMakeLists.txt` `project(gridcalc VERSION 0.2.0)` → `0.3.0`.
- Update `test/version_test.cpp` to assert `"0.3.0"` (rename of the existing `VersionTest.ReturnsZeroTwoZero` to `ReturnsZeroThreeZero`).
- `CHANGELOG.md`: new `## 0.3.0 — Phase 3: domain integration` entry summarizing the additions.
- `specs/STATUS.md`: bump "Last updated"; move Phase 3 row to Done; update "Where the project stands" + "Next action" (which becomes Phase 4 — Functional evaluation, callable API).

**Commit message:** `chore: bump version to 0.3.0 and refresh STATUS for Phase 3`
