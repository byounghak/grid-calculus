# Plan: Fix silent stencil aliasing on small periodic axes

## Group 1 — precondition helper + apply to `diff::` operators

- Add `include/gridcalc/diff/detail/PreconditionAxisExtent.hpp` with the
  free function
  `void requireAxisExtent(const char* axis, int N, int radius)` that
  throws `std::invalid_argument` with a precise message when
  `N < 2 * radius + 1`. Header carries `\since 0.14.1` Doxygen.
- Apply the precondition at the entry of every `diff::` operator that
  reads off-axis samples:
  - [include/gridcalc/diff/Laplacian.hpp](../../include/gridcalc/diff/Laplacian.hpp) — three axes always.
  - [include/gridcalc/diff/Gradient.hpp](../../include/gridcalc/diff/Gradient.hpp) — scalar overload (three axes), Phase 13 vector overload (three axes).
  - [include/gridcalc/diff/Divergence.hpp](../../include/gridcalc/diff/Divergence.hpp) — scalar-output overload (three axes), Phase 13 tensor-output overload (three axes).
  - [include/gridcalc/diff/Biharmonic.hpp](../../include/gridcalc/diff/Biharmonic.hpp) — three axes always.
  - [include/gridcalc/diff/detail/MultiIndexDerivative.hpp](../../include/gridcalc/diff/detail/MultiIndexDerivative.hpp) — only axes with the corresponding per-axis derivative count `> 0`. Single site for all 31 d-prefix mixed partials.
- Update each touched operator's Doxygen block with a `\throws std::invalid_argument` line referencing the precondition.

**Commit message:** `fix(diff): reject grids smaller than the stencil radius`

## Group 2 — apply precondition to `fvm::cellLaplacian`

- Add the same precondition to `fvm::cellLaplacian` in [include/gridcalc/fvm/CellLaplacian.hpp](../../include/gridcalc/fvm/CellLaplacian.hpp). Radius is `1` by construction at Phase 14, so this is currently a no-op for any `N ≥ 3`; keeps the operator-entry contract uniform across `diff::` and `fvm::`.
- Update the Doxygen `\throws` line on `fvm::cellLaplacian` accordingly.
- Update `fvm::diffuse` (the heterogeneous-D solver overload) Doxygen to note that the precondition propagates from `cellLaplacian`.

**Commit message:** `fix(fvm): forward stencil-radius precondition for defense-in-depth`

## Group 3 — tests

- New `test/stencil_axis_extent_test.cpp` covering:
  - **Throw tests** at `N = 2 * radius` (largest aliasing case) for each radius tier (1, 2, 3) across every affected operator: scalar Laplacian, scalar/vector gradient, scalar/tensor divergence, biharmonic, third- and fourth-derivative mixed partials at both orders, and FVM `cellLaplacian`.
  - **Per-axis label sweep** verifying the helper reports the right axis (x, y, z) in its message — three tests, one per axis, each tripping the check on a different axis with the others sized normally.
  - **Smallest-valid-N boundary tests** at `N = 2 * radius + 1` for every Phase 7 / Phase 8 / Phase 13 / Phase 14 operator. Each asserts no-throw and that the output is finite (no NaN / Inf). Rationale: extending the FD-FFT cross-check fixture's slope sweep to cover small-N is unsound — at `N = 5` the grid spacing `h = 2π/5 ≈ 1.257` is large enough that the truncation-error analysis underlying the slope test breaks down, and the slope-band assertion `[Order - 0.5, Order + 0.5]` would be polluted by the small-N point. The no-throw + finiteness gate at the boundary catches catastrophic alias-induced failures without making claims about truncation order at degenerate resolutions.
  - **Inactive-axis carve-out** verifying that an axis with per-axis derivative count `0` does not trip the precondition even at `N = 1`, through the `mixedDerivative<3, 0, 0, 2>` path on a grid with `Ny = Nz = 1`.

**Commit message:** `test(diff,fvm): cover stencil-radius precondition + smallest-valid-N boundary`

## Group 4 — version bump + CHANGELOG

- Bump `CMakeLists.txt` `VERSION` from `0.14.0` to `0.14.1`.
- Add a new `## 0.14.1 — Fix: silent stencil aliasing on small periodic axes` block at the top of `CHANGELOG.md` (above the `0.14.0` block) with `### Fixed` describing the bug + the precondition + the affected operator list; `### Tests` describing the new boundary coverage; `### Internal` referencing the new helper.
- Update `specs/STATUS.md`: add a "Phase 14 follow-up: PR #19" entry under the Phase 14 line, mirroring the Phase-10 follow-up format used for PR #14 (the Doxygen path-strip patch). Bump the **Last updated** date.

**Commit message:** `chore: release 0.14.1 — silent stencil aliasing fix`
