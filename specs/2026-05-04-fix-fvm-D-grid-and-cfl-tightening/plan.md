# Plan: Phase 14 hardening — D-grid mismatch + CFL tightening

## Group 1 — `fvm::cellLaplacian` preconditions (Grid-equality + D > 0)

- Edit [include/gridcalc/fvm/CellLaplacian.hpp](../../include/gridcalc/fvm/CellLaplacian.hpp):
  - Add `#include <gridcalc/solve/detail/RhsGridCheck.hpp>`.
  - At the entry of `cellLaplacian(psi, D)`:
    - `solve::detail::requireSameGrid(psi, D, "fvm::cellLaplacian: D")`.
    - O(N) scan over `D.data()`, throw `std::invalid_argument` on first cell with `D[idx] <= 0` (named flat index in the message).
  - Update Doxygen `\throws` line to enumerate the new preconditions.
  - `\since` tag: append `; 0.14.5 (D-grid + D > 0 preconditions)`.

**Commit message:** `fix(fvm): D-grid equality + D > 0 precondition at cellLaplacian entry`

## Group 2 — `solve::diffuse` updates (Grid-equality + drop redundant D > 0; CFL tightening)

- Edit [include/gridcalc/solve/Diffusion.hpp](../../include/gridcalc/solve/Diffusion.hpp):
  - Add `#include <gridcalc/solve/detail/RhsGridCheck.hpp>` if not already present.
  - At the entry of `diffuse(psi, D, dt, n_steps, tag)` (the heterogeneous-D overload):
    - `detail::requireSameGrid(psi, D, "solve::diffuse: D")`.
  - Drop the existing pre-loop `D > 0` scan (lines 163-176). Replace with a single O(N) pass computing `D_max_face = max over face pairs of fvm::detail::harmonicMean(D_a, D_b)`.
    - For each axis `a ∈ {x, y, z}` and each cell `(i, j, k)`, pair with `(i+1, j, k)` (or `(i, j+1, k)` / `(i, j, k+1)`) under periodic wrap.
    - Each face pair counted once across the grid (iterating only `+a` faces).
  - Pass `D_max_face` instead of `D_max` to `detail::checkDiffusionCFL<...>`.
  - Update Doxygen on the heterogeneous-D `solve::diffuse`:
    - `\throws` line: D-grid mismatch (new), D > 0 (now propagated from `cellLaplacian` rather than enforced here), CFL violation, etc.
    - Update the "CFL bound" passage to describe `D_max_face` as the spatial maximum of harmonic-mean face diffusivity (strictly ≤ `D_max`).

**Commit message:** `fix(solve): tighten heterogeneous-D CFL via max harmonic-mean face-D + D-grid check`

## Group 3 — tests

- New `test/fvm_grid_mismatch_test.cpp`:
  - **Grid-mismatch tests at `fvm::cellLaplacian`**: three flavors (different total cell count, same total / different per-axis extents, same shape / different cell size). Each asserts `EXPECT_THROW(..., std::invalid_argument)`.
  - **D > 0 throw test at `fvm::cellLaplacian`**: a D field with one zero / negative cell triggers a throw with a flat-index message.
  - **Grid-mismatch tests at `solve::diffuse(psi, D, ...)`**: same three flavors; error context names "solve::diffuse" so the user can tell which entry point caught the bug.
  - **D > 0 throw via `solve::diffuse`**: a D field with a non-positive cell now throws via `cellLaplacian` on the first stage; verify the error message format is sensible (no behavior regression for the user).
  - **Happy-path smoke**: a correctly-sized D field on `psi.getGrid()` integrates one step without throwing.
- New `test/fvm_cfl_tightening_test.cpp` (or append to existing `test/fvm_diffusion_acceptance_test.cpp`):
  - **High-contrast D**: one cell at `D = 100`, all neighbors at `D = 1`. Compute `D_max = 100`, `D_max_face ≈ 1.98`. Pick a `dt` such that `dt * 100 * Σ(1/h²) > 0.5` (would have thrown under the old bound) but `dt * 1.98 * Σ(1/h²) < 0.5` (passes under the new bound). Assert `solve::diffuse(...)` does **not** throw.
  - **Uniform D regression**: a uniform `D = D_const` field; the new `D_max_face` equals `D_const`. CFL behavior is byte-identical to the old `D_max` bound.
  - **Existing tests unchanged**: `test/heterogeneous_diffusion_test.cpp` and `test/fvm_diffusion_acceptance_test.cpp` use smooth varying D where `D_max_face ≈ D_max`; their dt choices stay valid under the new bound. Verify by running the full suite.

**Commit message:** `test(fvm): cover D-grid mismatch, D > 0 contract, and CFL tightening on high-contrast D`

## Group 4 — version bump + CHANGELOG + STATUS

- Bump `CMakeLists.txt` `VERSION` from `0.14.4` to `0.14.5`.
- Update `test/version_test.cpp` to assert `"0.14.5"`.
- Add a `## 0.14.5 — Fix: D-grid mismatch + CFL tightening on heterogeneous-D diffusion` block to `CHANGELOG.md` (above the `0.14.4` block) with `### Fixed`, `### Performance`, `### Tests` sections.
- Update `specs/STATUS.md`: add the merge-commit hash of PR #22 to the prior `0.14.4` entry; add a new "Phase 14 follow-up #5" entry; bump the **Last updated** line.

**Commit message:** `chore: release 0.14.5 -- D-grid + CFL tightening`
