# Plan: Phase 11 — Functional with up-to-4th-order gradients

## Group 1 — Extend `func::evaluate` to a 4-argument arity

- Edit `include/gridcalc/func/Functional.hpp`:
  - Add `#include <gridcalc/diff/Biharmonic.hpp>` to pull in Phase 8's biharmonic.
  - Insert a new `if constexpr` branch at the **top** of the dispatch ladder (so the longest-arity match wins under the existing tie-breaking rule):
    ```cpp
    if constexpr (std::is_invocable_r_v<double, F&,
                                        double, const core::Vec3d&,
                                        double, double>) {
      const auto grad   = diff::gradient(psi);
      const auto lap    = diff::laplacian(psi);
      const auto biharm = diff::biharmonic(psi);
      // pointwise loop over (Nx, Ny, Nz) materializing integrand(i,j,k) =
      //   f(psi(i,j,k), grad(i,j,k), lap(i,j,k), biharm(i,j,k))
    } else if constexpr (...) { /* existing 3-arg branch */ }
    ```
  - Update the function-level Doxygen block to document the fourth supported arity `f(double psi, const Vec3d& grad_psi, double lap_psi, double biharm_psi)` — `biharm_psi` is $\nabla^4\psi$ at the grid point as returned by `diff::biharmonic`.
  - Update the `static_assert` error message in the trailing `else` branch to enumerate all four supported arities (was three).
  - Adjust the `\since` tag to `0.4.0 (function); 0.11.0 (4-arg arity)`. Update the `\file` block similarly.
  - Existing 1-, 2-, and 3-arg call sites remain source-compatible (the new branch is appended at the head of the ladder, not replacing any existing one).

**Commit message:** `feat(func): extend evaluate() to 4-arg arity (psi, grad, lap, biharm)`

## Group 2 — PFC tests (hand-computed reference + convergence + spectral cross-check)

- Add `test/pfc_functional_test.cpp`:
  - Helper: `pfcDensity(double psi, double lap_psi, double biharm_psi, double q0, double eps)` returning the PFC integrand value, used by both the production-path test and the spectral cross-check (the integrand formula appears once in test code, not duplicated).
  - **`PfcFunctionalTest.HandComputedAtN64`** — `ψ = sin(x) sin(y) sin(z)` on `[0, 2π]³`, `q0 = 1.0`, `eps = 0.1`. Calls `func::evaluate(psi, pfc_lambda)` with `pfc_lambda = [&](double p, const Vec3d& g, double lp, double bp) { return pfcDensity(p, lp, bp, q0, eps); }`. Asserts `|F - F_ref| / |F_ref| < 2e-2` where `F_ref = 2.05546875 * π^3`. (`grad_psi` enters the lambda's signature so the 4-arg branch is selected, but is not used in the integrand — PFC has no gradient term.)
  - **`PfcFunctionalTest.SecondOrderConvergence`** — sweep `N ∈ {16, 32, 64}`, log-log slope of the relative error in `[1.7, 2.3]`, halving-`h` ratio in `[3.5, 4.5]`. Confirms the FD path converges at the documented `O(h²)` rate.
  - **`PfcFunctionalTest.SpectralCrossCheck`** — at `N = 32`, materialize `spectral::laplacian(psi)` and `spectral::biharmonic(psi)`, assemble the PFC integrand by hand into a `Field<double>`, reduce with `func::integrate`. Assert `|F_spectral - F_ref| / |F_ref| < 1e-12`. Confirms the analytical reference via an independent numerical path. Skipped (with `GTEST_SKIP`) when built with `-DGRIDCALC_USE_FFT=OFF`.
  - **`PfcFunctionalTest.PsiOnlyArityStillWorks`** — single-line regression test: `func::evaluate(psi, [](double p) { return p*p; })` returns `≈ π^3` to spectral precision. Verifies the new 4-arg branch did not accidentally shadow the 1-arg branch.
  - **`PfcFunctionalTest.GinzburgLandauStillWorks`** — single-line regression test reproducing the Phase 4 GL hand-computed value to confirm the 3-arg branch still dispatches correctly.
- Wire `pfc_functional_test.cpp` into `test/CMakeLists.txt`.
- Verify no other test file regresses.

**Commit message:** `test(func): add PFC hand-computed + convergence + spectral cross-check`

## Group 3 — Phase 11 doc chapters (User Guide + Developer Note)

Per Phase 10's "from Phase 11 onward, doc deliverables go straight to the LaTeX sources" decision, the doc deliverable lands directly under `docs/{user-guide,developer-note}/chapters/` — no markdown placeholders.

- **User Guide** (`docs/user-guide/chapters/`):
  - Add `11-higher-order-functional.tex` covering: the new 4-argument callable arity, when each argument is materialized, the PFC worked example walked through end-to-end (parameter choices → integrand structure → analytical $F_{\text{ref}}$ derivation → calling code → recovered numerical value at `N = 64` and the convergence sweep). This is the chapter the roadmap implies by "Worked example: phase-field-crystal free energy."
  - Renumber the existing placeholder chapters: `11-cahn-hilliard.tex` → `12-cahn-hilliard.tex`, `12-diamond-lattice.tex` → `13-diamond-lattice.tex`. Update `\input{}` lines in `docs/user-guide/main.tex`. The chapter-label anchors (`\label{chap:cahn-hilliard}`, `\label{chap:diamond-lattice}`) are unchanged.
- **Developer Note** (`docs/developer-note/chapters/`):
  - Add `10-higher-order-functional.tex` with the standard five-section structure (Theory · Math derivation · Algorithm · Design decisions · References):
    1. **Theory.** PFC free energy as a continuous functional, the role of the $(q_0^2 + \nabla^2)^2$ kernel in selecting a preferred wavelength (Swift–Hohenberg lineage), why a 4th-order *local* density requires materialized $\nabla^4\psi$ rather than $|\nabla(\nabla\psi)|^2$ (the two are not equal at finite `h`).
    2. **Math derivation.** Expansion of $(q_0^2+\nabla^2)^2 = q_0^4 + 2q_0^2\nabla^2 + \nabla^4$; closed-form $F_{\text{ref}} = 2.05546875\,\pi^3$ for the test case; leading-order discrete error scaling argument $|F - F_h| \sim h^2$ from the biharmonic stencil.
    3. **Algorithm.** SFINAE arity-dispatch ladder extended at the head; eager materialization of `gradient`, `laplacian`, `biharmonic`; reuse of Phase 8's fused 6-term direct biharmonic stencil (not `laplacian∘laplacian`); algorithmic complexity (one extra grid pass for biharmonic, one extra `Field<double>` allocation for its result).
    4. **Design decisions.** Cite `requirements.md` for: contracted-scalar-only argument shape; SFINAE-arity dispatch over variadic context object; deferral of contraction ETs to Phase 14; spectral cross-check rather than FD–FFT-fixture extension.
    5. **References.** Elder & Grant (2004), *Modeling elastic and plastic deformations in nonequilibrium processing using phase field crystals*, Phys. Rev. E 70, 051605 — DOI link for the canonical PFC functional. Provatas & Elder, *Phase-Field Methods in Materials Science and Engineering* (textbook, chapter on PFC). Boyd, *Chebyshev and Fourier Spectral Methods*, 2nd ed. (textbook, on spectral exactness for band-limited inputs — used in the spectral cross-check).
  - Update `docs/developer-note/main.tex` to `\input{chapters/10-higher-order-functional.tex}` after the Phase 9 chapter, before `\appendix`.
- Both new chapters carry chapter labels (`\label{chap:higher-order-functional}` etc.) so other chapters can `\Cref{}` them later.

**Commit message:** `docs: add Phase 11 User Guide and Developer Note chapters`

## Group 4 — Bookkeeping

- Bump `CMakeLists.txt`: `project(gridcalc VERSION 0.10.0 ...)` → `0.11.0`.
- Update `test/version_test.cpp` accordingly (the existing version-string check).
- Bump the `Version \texttt{0.10.0}` strings on the title pages of both `docs/user-guide/main.tex` and `docs/developer-note/main.tex` to `0.11.0`.
- `CHANGELOG.md`: new `## 0.11.0 — Phase 11: Functional with up-to-4th-order gradients` entry under the standard Keep-a-Changelog headings (`Added` / `Changed`).
- `specs/STATUS.md`: bump **Last updated** to today; move Phase 11 row to Done; update "Where the project stands" with a Phase-11 paragraph + PR link placeholder; refresh "Next action" to Phase 12 (Cahn–Hilliard end-to-end demo); add the contracted-scalar-only and ET-deferral decisions to "Decisions worth knowing" if they are not yet captured there.

**Commit message:** `chore: bump version to 0.11.0 and refresh STATUS for Phase 11`
