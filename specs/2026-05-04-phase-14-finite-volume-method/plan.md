# Plan: Phase 14 — Finite volume method (FVM)

Eight commit-sized groups. Group 7 was added mid-phase per a user request to ship a textbook FVM demo. Build and tests stay green between every commit.

## Group 1 — Spec files

- `specs/2026-05-04-phase-14-finite-volume-method/{requirements, plan, validation}.md`.

**Commit message:** `docs: phase 14 spec files (finite volume method)`

## Group 2 — Roadmap insertion + mechanical renumber

- `specs/roadmap.md` — insert new Phase 14 entry between current Phase 13 and current (about-to-be-renamed) Phase 14; renumber `Phase N → Phase N+1` for `N ∈ [14, 22]`. Update internal cross-references (e.g., the Phase 21 / Phase 22 mentions in the Phase 0 acceptance bullet, the Phase 16 / Phase 22 mentions in the Phase 10 goal, the Risks-to-watch section).
- `specs/STATUS.md` — "Next action" rewritten to point to new Phase 14 (FVM); "Decisions worth knowing" cross-references renumbered via a Python regex pass constrained to the relevant section (`Phase N → Phase N+1` for `N ∈ [14, 22]`, descending); "Phase progress" table updated.
- `docs/user-guide/chapters/13-vector-tensor-fields.tex` — `Phase 14` → `Phase 15` (the deferred-tensor-functional-eval reference).
- `docs/developer-note/chapters/12-vector-tensor-fields.tex` — `Phase 14` → `Phase 15` (multiple references).
- `docs/user-guide/chapters/15-diamond-lattice.tex` — `Phase 16` → `Phase 17`.
- `test/integrate_test.cpp` — comment `Phase 20` → `Phase 21`.
- Sanity-check by grep that no other live references to "Phase 1[4-9]\|Phase 2[0-2]" remain stale.

CHANGELOG entries for 0.11.0 / 0.12.0 / 0.13.0 are **not** renumbered (historical records, git-immutable). The new 0.14.0 entry in Group 8 explicitly calls out the renumber.

**Commit message:** `docs: insert Phase 14 (FVM); renumber existing phases 14-22 to 15-23`

## Group 3 — `gridcalc::fvm::cellLaplacian` (constant-D agreement gate)

- New file `include/gridcalc/fvm/CellLaplacian.hpp` — exposes `Field<double> fvm::cellLaplacian(const Field<double>& psi, const Field<double>& D)` returning `∇·(D ∇ψ)` via cell-flux discretization with face-centered harmonic-mean D-averaging. Doxygen block with `\since 0.14.0`, the explicit face-flux formula, the harmonic-mean choice rationale, and the contract that `D > 0` everywhere. The detail helper `harmonicMean(a, b)` lives in `gridcalc::fvm::detail`.
- `test/fvm_cell_laplacian_test.cpp` — first three tests:
  - `AgreesWithFdLaplacianAtConstantUnitD` — `fvm::cellLaplacian(psi, D=1.0_field)` matches `diff::laplacian(psi)` to round-off (max-abs < 1e-12) on a 16³ trig manufactured solution.
  - `AgreesWithFdLaplacianAtConstantD` — same with `D = c_field` for `c = 0.7`; result equals `c · diff::laplacian(psi)`.
  - `ConstantPsiGivesZeroRhs` — `psi = const` ⇒ RHS is zero (FVM is conservative).

**Commit message:** `feat(fvm): add cellLaplacian heterogeneous-D operator (constant-D agreement gate)`

## Group 4 — Heterogeneous-D unit tests

- `test/fvm_cell_laplacian_test.cpp` — extends with three more tests:
  - `HarmonicMeanFaceCoefficient` — 4-cell strip with `D = [1, 4, 4, 1]` and `psi = [0, 1, 1, 0]`; verifies the face flux at the D-jump cells matches the harmonic-mean prediction analytically (e.g., `hmean(1, 4) = 1.6`, so cell-1 RHS = -1.6).
  - `AnalyticalSolutionConvergence_Order2` — manufactured solution `D(x) = 1 + 0.5 cos(x)`, `ψ(x) = sin(x)`. Closed-form `∇·(D ∇ψ) = -sin(x) - 0.5 sin(2x)`. Convergence sweep on `N ∈ {16, 32, 64}` asserts slope ∈ `[1.7, 2.3]`.
  - `MassConservation_PeriodicSum` — for arbitrary `D > 0` and `ψ`, `Σ fvm::cellLaplacian(ψ, D) = 0` to round-off (telescoping fluxes on a periodic domain).

**Commit message:** `test(fvm): heterogeneous-D convergence + harmonic-mean + mass conservation`

## Group 5 — `solve::diffuse` heterogeneous-D overload

- `include/gridcalc/solve/Diffusion.hpp` — extend with the new overload `template <class Tag> void diffuse(Field<double>& psi, const Field<double>& D, double dt, int n_steps, Tag tag)`. Pre-loop pass over `D` validates the harmonic-mean contract (`D > 0` everywhere) and computes `D_max` for the CFL bound in the same sweep. CFL bound generalizes to `D_max · dt · sum_a (1/h_a²) ≤ Tag::diffusionCFLLimit`. RHS callable wraps `fvm::cellLaplacian(psi, D)`. The `D` factor lives inside the operator (FVM returns `∇·(D ∇ψ)` directly), so the integrator step is just `dt` — no `D · dt` folding like the constant-D path. `\since 0.14.0`.
- `test/diffusion_test.cpp` — extend with four new tests:
  - `HeterogeneousD_AgreesWithConstantDOnUniformField` — feeding the heterogeneous-D path a uniform `D = c` field reproduces the constant-D path's trajectory to round-off.
  - `HeterogeneousD_RejectsCFLViolation` — overshooting the generalized bound throws `std::invalid_argument`.
  - `HeterogeneousD_RejectsNonpositiveD` — `D ≤ 0` at any cell throws (with a flat-index pointer in the message). Both negative and zero values are rejected.
  - `HeterogeneousD_MassConservedOverTimeStepping` — Gaussian IC diffused through a heterogeneous `D = 0.05 + 0.04 sin(x)` field for 100 explicit-Euler steps; total mass conserved to round-off.

**Commit message:** `feat(solve): heterogeneous-D solve::diffuse overload reusing existing tags`

## Group 6 — End-to-end acceptance test

- `test/fvm_diffusion_acceptance_test.cpp` — new file. Two tests discharge the roadmap acceptance bar:
  - `RecoversAnalyticalEigenfunctionAtOrder2` — runs the heterogeneous-D `solve::diffuse` code path with a uniform `D = 0.1` field on the well-known `sin(x) sin(y) sin(z)` eigenfunction. The trajectory must match the closed-form `exp(-3 D t) · sin(x) sin(y) sin(z)` at order 2 in `h` on `N ∈ {16, 24, 32}`; slope ∈ `[1.6, 2.4]`. Exercises every line of the new path (cellLaplacian with harmonic-mean averaging, the heterogeneous-D `solve::diffuse` overload, the CFL check, the D-positivity validation, the integrator loop) on a problem with a known analytical answer.
  - `HeterogeneousDQualitativeProperties` — runs with a genuinely varying `D(x) = 0.10 + 0.05 cos(x)` (range [0.05, 0.15]) and a Gaussian peak IC; runs 200 RK4 steps in 4 chunks and asserts three properties that hold for any `D > 0`: mass conservation, max|ψ| monotonically non-increasing across chunks, and `ψ ≥ 0` throughout (positivity preservation).

  **Why two tests, not one.** The naive manufactured solution `ψ(x, y, z, t) = exp(-D(x) t) sin(z)` only satisfies the PDE at `t = 0`; once ψ acquires x-dependence from the position-dependent decay rate, the x-axis flux is non-zero and the analytical reference drifts. A simple closed-form solution to the genuinely heterogeneous-D PDE on a periodic domain is not readily available, so the two-test combination above sidesteps the issue with the eigenfunction (uniform-D path, analytical reference) plus the qualitative gates (heterogeneous-D path, no analytical reference but tight property assertions). The trap is documented in the file's top-of-file comment.

**Commit message:** `test(fvm): heterogeneous-D end-to-end acceptance — 2nd-order convergence + qualitative properties`

## Group 7 — Textbook FVM demo (`examples/heterogeneous_diffusion.cpp`) — added 2026-05-04

Mid-phase user request: ship a runnable textbook nontrivial FVM example mirroring the Phase 12 CH demo's pattern.

- `examples/common/heterogeneous_diffusion.hpp` — header-only helper exposing `makeGaussianIc(grid, sigma, center)`, `makeHeterogeneousD(grid, D0, eps, axis)` (smooth `D = D0 + eps * cos(2π·k_axis/L)` on a chosen axis), and energy diagnostics (`computeMass(psi)`, `computeL2Squared(psi)`, `computeMin(psi)`, `computeMax(psi)`).
- `examples/heterogeneous_diffusion.cpp` — CLI demo. CLI options: `--n-x`, `--dt`, `--n-steps`, `--snapshot-every`, `--d-mean`, `--d-amplitude`, `--gaussian-sigma`, `--out-dir`, `--integrator {euler|rk4}`. Produces `psi_NNNNNN.vtk` snapshots via the Phase 12 `examples/common/vtk_writer.hpp`. Per-snapshot diagnostics printed to stdout (mass, L²-norm, max, min).
- `test/heterogeneous_diffusion_test.cpp` — small CI gate (~20–30 s budget) that runs a stripped-down version of the demo at `N=24`, `n_steps=200` and asserts mass conservation + L²-norm decrease + positivity preservation. Mirrors the Phase 12 `cahn_hilliard_test.cpp` static-cached-snapshot pattern so the simulation runs once per process.

**Commit message:** `feat(examples): textbook heterogeneous-D diffusion demo + CI gate`

## Group 8 — Doc deliverables and bookkeeping

- `docs/user-guide/chapters/14-finite-volume-method.tex` — full chapter (slot 14 is currently empty; chapter 15 is the diamond-lattice placeholder for the renamed Phase 16). Sections: motivation (when to use FVM vs. FD); the cell-flux discretization with the explicit face-flux formula; harmonic-mean face-D averaging derivation (steady-state-flux match across a discontinuity); the heterogeneous `solve::diffuse` API; the textbook demo walkthrough (Group 7) with three committed PNG snapshots; what this does *not* do.
- `docs/developer-note/chapters/13-finite-volume-method.tex` — five-section structure (Theory · Math derivation · Algorithm · Design decisions · References). Theory: conservation laws and the flux-form PDE; Bray-style framing. Math derivation: the cell-volume integration, face-flux approximation, harmonic-mean derivation, order-2 truncation analysis, generalized CFL. Algorithm: file layout, the `cellLaplacian` loop, the heterogeneous-D `solve::diffuse` overload, the demo's IC + diagnostics helpers. Design decisions: the four AskUserQuestion answers + the manufactured-solution trap + the no-FD-replacement rationale + the demo addition. References: LeVeque 2007 (FVM section), Patankar 1980, Hundsdorfer & Verwer 2003, Eigen project docs.
- `docs/user-guide/figures/heterogeneous-diffusion/{early,mid,late}.png` — three snapshot PNGs from a 64³ release-build run via a new `scripts/render_diffusion_snapshots.py` (reuses the Phase 12 matplotlib "Agg" pattern).
- `docs/user-guide/main.tex` and `docs/developer-note/main.tex` — add `\input{}` lines.
- Bump `project(... VERSION 0.14.0)` in root `CMakeLists.txt`.
- Update `test/version_test.cpp` to assert `0.14.0`.
- New `## 0.14.0 — Phase 14` block in `CHANGELOG.md` (calls out the roadmap renumber explicitly).
- Refresh `specs/STATUS.md`: Phase 14 row → Done; "Last updated" date bump; "Next action" rewritten for Phase 15 (the renumbered functional-evaluation-for-vector/tensor-data); new decisions added.

**Commit message:** `chore: bump version to 0.14.0 and refresh STATUS for Phase 14`

---

After Group 8: `git push -u origin 2026-05-04-phase-14-finite-volume-method`, open PR, watch CI per `specs/CLAUDE.md` step 6.
