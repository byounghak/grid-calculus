# Plan: Phase 14 — Finite volume method (FVM)

Seven commit-sized groups. Build and tests stay green between every commit.

## Group 1 — Spec files

- `specs/2026-05-04-phase-14-finite-volume-method/{requirements, plan, validation}.md`.

**Commit message:** `docs: phase 14 spec files (finite volume method)`

## Group 2 — Roadmap insertion + mechanical renumber

- `specs/roadmap.md` — insert new Phase 14 entry between current Phase 13 and current (about-to-be-renamed) Phase 14; renumber `Phase N → Phase N+1` for `N ∈ [14, 22]`. Update internal cross-references (e.g., the Phase 21 / Phase 22 mentions in the Phase 0 acceptance bullet).
- `specs/STATUS.md` — "Next action" rewritten to point to new Phase 14 (FVM); "Decisions worth knowing" cross-references renumbered (`Phase 14` → `Phase 15`, `Phase 18` → `Phase 19`, `Phase 19` → `Phase 20`, `Phase 20` → `Phase 21`, `Phase 21` → `Phase 22`); "Phase progress" table updated.
- `docs/user-guide/chapters/13-vector-tensor-fields.tex` — `Phase 14` → `Phase 15` (the deferred-tensor-functional-eval reference).
- `docs/developer-note/chapters/12-vector-tensor-fields.tex` — `Phase 14` → `Phase 15` (multiple references).
- `test/integrate_test.cpp` — comment `Phase 20` → `Phase 21`.
- Sanity-check by grep that no other live references to "Phase 1[4-9]\|Phase 2[0-2]" remain stale.

CHANGELOG entries for 0.11.0 / 0.12.0 / 0.13.0 are **not** renumbered (historical records, git-immutable). The new 0.14.0 entry in Group 7 explicitly calls out the renumber.

**Commit message:** `docs: insert Phase 14 (FVM); renumber existing phases 14-22 to 15-23`

## Group 3 — `gridcalc::fvm::cellLaplacian` (constant-D agreement gate)

- New file `include/gridcalc/fvm/CellLaplacian.hpp` — exposes `Field<double> fvm::cellLaplacian(const Field<double>& psi, const Field<double>& D)` returning `∇·(D ∇ψ)` via cell-flux discretization with face-centered harmonic-mean D-averaging. Doxygen block with `\since 0.14.0`, the explicit face-flux formula, the harmonic-mean choice rationale, and the contract that `D > 0` everywhere.
- `test/fvm_cell_laplacian_test.cpp` — first three tests:
  - `AgreesWithFdLaplacianAtConstantUnitD` — `fvm::cellLaplacian(psi, D=1.0_field)` matches `diff::laplacian(psi)` to round-off on a periodic trig manufactured solution.
  - `AgreesWithFdLaplacianAtConstantD` — same but `D = c_field` for `c = 0.7`; result should equal `c · diff::laplacian(psi)`.
  - `ConstantPsiGivesZeroRhs` — `psi = const` ⇒ RHS is zero (FVM is conservative; constant ψ has zero flux).

**Commit message:** `feat(fvm): add cellLaplacian heterogeneous-D operator (constant-D agreement gate)`

## Group 4 — Heterogeneous-D tests

- `test/fvm_cell_laplacian_test.cpp` — extends with three more tests:
  - `HarmonicMeanFaceCoefficient` — constructs a 2-cell test where `D` jumps from `D_low` to `D_high` across one face; verifies the face flux matches the harmonic-mean prediction analytically.
  - `AnalyticalSolutionConvergence_Order2` — for a manufactured `D(x) = 1 + 0.5 cos(x)` and `ψ(x) = sin(x)`, the closed-form `∇·(D ∇ψ)` is `−sin(x)·(1 + 0.5 cos(x)) + 0.5 sin(x)·cos(x)·(−1)`. Convergence sweep on `N ∈ {16, 32, 64}` asserts slope `∈ [1.7, 2.3]`.
  - `MassConservation_PeriodicSum` — for arbitrary `D > 0` and `ψ`, `Σ fvm::cellLaplacian(ψ, D) = 0` to round-off (telescoping fluxes on a periodic domain).

**Commit message:** `test(fvm): heterogeneous-D convergence + harmonic-mean + mass conservation`

## Group 5 — `solve::diffuse` heterogeneous-D overload

- `include/gridcalc/solve/Diffusion.hpp` — extend with the new overload `template <class Tag> void diffuse(Field<double>& psi, const Field<double>& D, double dt, int n_steps, Tag tag)`. Body: validates the CFL bound `D_max · dt · Σ_a (1/h_a²) ≤ Tag::diffusionCFLLimit`, then calls `solve::integrate(psi, [&D](const auto& y){ return fvm::cellLaplacian(y, D); }, dt, n_steps, tag)`. Constant-D Phase 5 / Phase 6 overload unchanged. `\since 0.14.0`.
- `test/diffusion_test.cpp` — extend with three new tests:
  - `HeterogeneousD_AgreesWithConstantDOnUniformField` — running both the constant-D and the heterogeneous-D paths with the same uniform `D` produces the same trajectory after `n_steps` to round-off.
  - `HeterogeneousD_RejectsCFLViolation` — passing a `dt` above the heterogeneous CFL bound throws `std::invalid_argument`.
  - `HeterogeneousD_MassConservedOverTimeStepping` — `Σ ψ` is conserved to round-off across many steps (the FVM property carrying through the integrator).

**Commit message:** `feat(solve): heterogeneous-D solve::diffuse overload reusing existing tags`

## Group 6 — Acceptance test (manufactured-solution end-to-end)

- `test/fvm_diffusion_acceptance_test.cpp` — new file. Discharges the roadmap acceptance bar: with `D(x) = 1 + 0.5 cos(x)`, `ψ_0(x) = sin(x)`, the closed-form 1D heterogeneous-diffusion equation has a known time-decay; verify the discrete trajectory matches the analytical solution at `t = 1.0` to within `O(h²)` tolerance on `N ∈ {32, 64}`. Uses RK4 to keep time-discretization error well below the spatial-discretization error so the spatial-convergence rate is observable.

**Commit message:** `test(fvm): heterogeneous-D end-to-end manufactured-solution acceptance`

## Group 7 — Doc deliverables and bookkeeping

- `docs/user-guide/chapters/14-finite-volume-method.tex` — full chapter (slot 14 is currently empty; chapter 15 is the diamond-lattice placeholder for the renamed Phase 16). Sections: motivation (when to use FVM vs. FD), the cell-flux discretization, harmonic-mean face-D averaging, the heterogeneous-`solve::diffuse` overload, the manufactured-solution acceptance test as a code example.
- `docs/developer-note/chapters/13-finite-volume-method.tex` — five-section structure (Theory · Math derivation · Algorithm · Design decisions · References). Theory: conservation laws and the flux-form PDE; Bray-class advection-diffusion textbook framing. Math derivation: the cell-volume integration, face-flux approximation, harmonic-mean derivation (steady-state flux match across a discontinuity), order-2 truncation analysis. Algorithm: file layout, the `cellLaplacian` loop, the CFL bound generalization. Design decisions: AskUserQuestion answers + the no-FD-replacement rationale. References: LeVeque 2007 (FVM section), Patankar 1980, Hundsdorfer & Verwer 2003 (ADR equation reference), Eigen project docs.
- `docs/user-guide/main.tex` — add `\input{chapters/14-finite-volume-method.tex}` between chapter 13 and chapter 15.
- `docs/developer-note/main.tex` — add `\input{chapters/13-finite-volume-method.tex}` after chapter 12.
- Bump `project(... VERSION 0.14.0)` in root `CMakeLists.txt`.
- Update `test/version_test.cpp` to assert `0.14.0`.
- New `## 0.14.0 — Phase 14` block in `CHANGELOG.md` (calls out the roadmap renumber explicitly so future readers can correlate "Phase 14" in earlier blocks with the renamed phase).
- Refresh `specs/STATUS.md`: Phase 14 row (FVM) → Done; "Last updated" date bump; "Next action" rewritten for Phase 15 (the renumbered functional-evaluation-for-vector/tensor-data); new decisions added.

**Commit message:** `chore: bump version to 0.14.0 and refresh STATUS for Phase 14`

---

After Group 7: `git push -u origin 2026-05-04-phase-14-finite-volume-method`, open PR, watch CI per `specs/CLAUDE.md` step 6.
