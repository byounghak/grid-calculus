# Plan: Phase 12 — Cahn–Hilliard example end-to-end

Seven commit-sized groups. The build and the test suite stay green between every commit.

## Group 1 — Spec files

- `specs/2026-05-03-phase-12-cahn-hilliard/requirements.md`
- `specs/2026-05-03-phase-12-cahn-hilliard/plan.md` (this file)
- `specs/2026-05-03-phase-12-cahn-hilliard/validation.md`

**Commit message:** `docs: phase 12 spec files (Cahn–Hilliard example)`

## Group 2 — Examples build wiring

- Add top-level option `GRIDCALC_BUILD_EXAMPLES` (default `OFF`) to root `CMakeLists.txt`, mirroring `GRIDCALC_BUILD_DOCS`.
- `examples/CMakeLists.txt` gating the example targets behind the option.
- Placeholder `examples/cahn_hilliard.cpp` with `int main() { return 0; }` so the wiring builds.
- Update both Ubuntu CI jobs in `.github/workflows/ci.yml` to pass `-DGRIDCALC_BUILD_EXAMPLES=ON` so the example compiles on every PR (not run; just type-checked).

**Commit message:** `build: add examples/ subtree and GRIDCALC_BUILD_EXAMPLES option`

## Group 3 — VTK writer helper

- `examples/common/vtk_writer.hpp` — single header, single function `writeVtkStructuredPoints(const Field<double>&, const Grid&, const std::string& path, const std::string& field_name)`. Writes legacy ASCII VTK 3.0 STRUCTURED_POINTS.
- `test/vtk_writer_test.cpp` — round-trip coverage:
  - `WriteThenParse_HeaderMatchesGrid` — write a small 4×3×2 field, parse the header back, check `DIMENSIONS` / `ORIGIN` / `SPACING` / `POINT_DATA` line counts.
  - `WriteThenParse_PayloadMatchesField` — parse the floats and assert exact equality with the source field's flat layout.
  - `WriteThenParse_FieldNameRoundTrips` — the `SCALARS <name> double 1` line uses the user-supplied name.

**Commit message:** `feat(examples): add legacy ASCII VTK STRUCTURED_POINTS writer`

## Group 4 — CH RHS + main

- `examples/common/cahn_hilliard.hpp` — exposes `Field<double> computeCahnHilliardRhs(const Field<double>& psi, double M, double kappa)`. One-liner body: `M * laplacian(-psi + psi³ - kappa * laplacian(psi))` realized as Phase-1 / Phase-2 calls with the appropriate pointwise materializations.
- `examples/cahn_hilliard.cpp` — full implementation. CLI parsing for `--n-x`, `--dt`, `--n-steps`, `--snapshot-every`, `--kappa`, `--mobility`, `--seed`, `--out-dir`. Uniform random IC ψ ~ U(-0.05, 0.05) with empirical mean subtracted. Time-stepping loop: in chunks of `--snapshot-every` steps, call `solve::integrate(psi, rhs, dt, snapshot_every, RK4{})` then write a VTK snapshot via `writeVtkStructuredPoints`. Print free energy + gradient energy + mean to `stdout` per snapshot.

**Commit message:** `feat(examples): Cahn–Hilliard coarsening demo with RK4 + VTK output`

## Group 5 — Acceptance test

- `test/cahn_hilliard_test.cpp` — five tests:
  - `RhsClosedFormOnSinusoid` — for ψ = sin(x)sin(y)sin(z) on a small grid, compare `computeCahnHilliardRhs` against the hand-computed analytic RHS within stencil-order error.
  - `LyapunovDecay` — run `N = 32`, `~3000` RK4 steps with random IC, take `K = 6` snapshots; assert F(t_{i+1}) ≤ F(t_i) + tol_F at every snapshot.
  - `MonotoneCoarsening` — same run, gradient-energy `∫|∇ψ|²` decreases between the post-transient snapshot pair (e.g. snapshot 4 to snapshot 5).
  - `MassConservation` — same run, `|⟨ψ⟩(t_final) − ⟨ψ⟩(0)| < 1e-10`.
  - `SeedReproducibility` — running the demo's IC generator with a fixed seed twice yields bit-identical fields.

**Commit message:** `test(examples): assert CH free-energy decay + gradient-energy monotone coarsening`

## Group 6 — Doc deliverables + figures

- `docs/user-guide/chapters/12-cahn-hilliard.tex` — replace placeholder with the full chapter (equation → code → results), per the User Guide style established in Phases 10–11. Sections: 1 Spinodal decomposition / motivation; 2 Free-energy form and variational derivative; 3 Discretization choices and RK4 stability budget; 4 Code walkthrough; 5 Running the demo + VTK output; 6 Visual results (3 PNGs); 7 Where to go next (Phase 19 implicit treatment, Phase 13/14 vector/tensor extension).
- `docs/developer-note/chapters/11-cahn-hilliard.tex` — five-section structure (Theory · Math derivation · Algorithm · Design decisions · References). Theory: Cahn 1958, Cahn–Hilliard 1958, the Lyapunov property dF/dt = -M ∫ |∇μ|² ≤ 0. Math derivation: the closed-form δF/δψ; the linear stability analysis σ(k) = M k²(1 - κk²); the explicit-RK4 stability radius applied to the worst-case eigenvalue M·κ·k_max⁴. Algorithm: file-paths + control-flow. Design decisions: the four AskUserQuestion answers + the rejected alternatives. References: at least Cahn–Hilliard 1958, Bray 1994, Provatas–Elder 2010, LeVeque 2007 (RK4 stability).
- Update `docs/user-guide/main.tex` and `docs/developer-note/main.tex` if needed (likely no change — placeholder chapter was already `\input`'d).
- `docs/user-guide/figures/cahn-hilliard/{early,mid,late}.png` — 2D z-midplane slices.
- `scripts/render_ch_snapshots.py` — matplotlib helper that reads the example's `.vtk` output and writes the three PNGs.

**Commit message:** `docs: User Guide ch.12 + Developer Note ch.11 for Cahn–Hilliard`

## Group 7 — Bookkeeping

- Bump `project(... VERSION 0.12.0)` in root `CMakeLists.txt`.
- `CHANGELOG.md` — new `## 0.12.0 — Phase 12` block following the existing format.
- `specs/STATUS.md` — Phase 12 row to "Done", date bump, "Next action" rewrites for Phase 13 (vector/tensor fields), new decisions added to "Decisions worth knowing", `MEMORY.md`-style PR list extended.
- `specs/roadmap.md` — no changes expected; reread to confirm Phase 12 entry still matches what we shipped.

**Commit message:** `chore: bump version to 0.12.0 and refresh STATUS for Phase 12`

---

After Group 7: `git push -u origin 2026-05-03-phase-12-cahn-hilliard`, open the PR, watch CI per `specs/CLAUDE.md` step 6.
