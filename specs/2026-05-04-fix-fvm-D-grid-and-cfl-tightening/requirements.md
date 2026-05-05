# Requirements: Phase 14 hardening — D-grid mismatch + CFL tightening for heterogeneous-D diffusion

## Bug report

Two issues in the Phase 14 (FVM / heterogeneous-D diffusion) surface:

### Issue 1: missing D-grid equality and D > 0 checks at the operator entry

[`solve::diffuse(psi, D, ...)`](../../include/gridcalc/solve/Diffusion.hpp#L148-L185)
validates `D > 0` (good) and computes `D_max` for the CFL bound, but
**does not check `D.getGrid() == psi.getGrid()`**. If `D` was sampled
on a different grid than `psi`:

- **Different extents, same total size**: `D(i, j, k)` returns a value
  at a different physical location than `psi` expects. Fluxes are
  wrong, no diagnostic.
- **Different cell size `h`**: `cellLaplacian` uses `psi`'s `h` for
  `inv_hx2` etc.; if `D` was sampled on a grid with different `h`,
  the spatial scale of `D` vs `psi` diverges silently.
- **Smaller `D`**: `cellLaplacian` indexes `D(i ± 1, ...)` for `psi`'s
  full extent; `D`'s own periodic wrap masks OOB but produces nonsense
  values.

[`fvm::cellLaplacian(psi, D)`](../../include/gridcalc/fvm/CellLaplacian.hpp#L75-L111)
itself has **no Grid-equality check** at the operator entry, and **no
`D > 0` validation** — that contract is enforced only by the
`solve::diffuse` wrapper. Direct callers of `cellLaplacian` (any user
code that wraps it as the RHS for their own time stepper) get UB /
NaN / Inf if they violate either contract. Same defensive-contract
gap as PR #20 closed for `solve::integrate` / `axpyFresh`.

### Issue 2: CFL bound is conservative for high-contrast `D` (perf cliff)

`solve::diffuse(psi, D, ...)` computes `D_max = max over cells of D`
and uses it for the von-Neumann CFL bound:
`D_max · dt · Σ_a (1/h_a²) ≤ Tag::diffusionCFLLimit`.

But the cell-flux operator uses **harmonic-mean face diffusivity**,
which is bounded above by `min(D_a, D_b)` (with equality at
`D_a = D_b`). For high-contrast `D` (e.g., one cell at `D = 100`
surrounded by cells at `D = 1`), the harmonic mean strongly suppresses
the maximum face-D: `harmonicMean(1, 100) = 200/101 ≈ 1.98`, not 100.
Using `D_max = 100` forces an unnecessarily tiny `dt` — a ~50× perf
cliff for high-contrast cases.

**Not a correctness bug** — strictly conservative — but a perf /
usability cliff that disappears if `D_max` is replaced with
`D_max_face = max over face pairs of harmonicMean(D_a, D_b)`. For
uniform `D`, `harmonicMean(D, D) = D`, so existing constant-D and
smooth-D users see byte-identical trajectories.

## Why dormant in this codebase (Issue 1)

- The `fvm::cellLaplacian` overload of `solve::diffuse` builds its
  RHS lambda as `[&D](const Field<double>& f) { return fvm::cellLaplacian(f, D); }`,
  capturing the user-passed `D`. The user is responsible for passing
  a `D` allocated on the same `Grid` as `psi`.
- The repo's own tests (`test/heterogeneous_diffusion_test.cpp`,
  `test/fvm_diffusion_acceptance_test.cpp`,
  `examples/heterogeneous_diffusion.cpp`) construct `D` on
  `psi.getGrid()` directly — no mismatch is ever exercised.
- Reachable by user code that holds a captured `Grid` for `D` (e.g.,
  loaded from a different source) or that composes `cellLaplacian`
  outside `solve::diffuse`.

## Decisions made for this feature

- **Failure mode: `throw std::invalid_argument`** at the operator
  entry. Same pattern as the prior three patches
  (`requireAxisExtent`, `requireSameGrid`, etc.).

- **Reuse the existing `solve::detail::requireSameGrid` helper.**
  Already in place from PR #20 (`0.14.2`). Cross-namespace include
  from `fvm::` to `solve::detail::` is consistent with the existing
  `fvm::cellLaplacian` → `diff::detail::requireAxisExtent` reuse.
  No relocation of the helper to `core::detail::` — keeps the patch
  surface small and matches the established pattern.

- **D > 0 validation lives at `fvm::cellLaplacian` entry only.**
  Direct callers of `cellLaplacian` get the loud failure. The
  current `solve::diffuse` pre-loop scan that combined `D > 0` with
  `D_max` computation is dropped — `cellLaplacian` now owns the
  positivity contract per-call. Total work per `solve::diffuse` call
  is unchanged: `solve::diffuse`'s pre-loop pass still runs once
  (now computing `D_max_face` only), and `cellLaplacian`'s per-call
  scan replaces what the pre-loop scan used to do — net cost is
  the same `O(N · n_steps)`. The difference is that the contract
  is now enforced at the operator (so direct callers benefit) and
  the error message context shifts from "solve::diffuse: pre-loop"
  to "fvm::cellLaplacian: D".

- **Grid-equality check at both entry points.** `solve::diffuse`
  validates up front so the user gets a clearer "solve::diffuse: D"
  context before the integrator starts; `cellLaplacian` validates
  per-call so direct callers benefit. The cellLaplacian-internal
  check is technically redundant when called from solve::diffuse,
  but it's `O(1)` (a single `Grid::operator==` comparison), so the
  overhead is negligible.

- **CFL bound: `D_max_face = max over face pairs of harmonicMean(D_a, D_b)`.**
  Computed in a single O(N) scan over `+x / +y / +z` face pairs (the
  `−x / −y / −z` faces are the same pairs viewed from the adjacent
  cell, so iterating only `+` faces is sufficient). Bit-identical to
  the current `D_max` bound for uniform `D` (since
  `harmonicMean(D, D) = D`). Strictly tighter for varying `D`.

- **No CFL helper signature change.** `detail::checkDiffusionCFL`
  still takes `(grid, double D_eff, dt)` — we just pass
  `D_max_face` instead of `D_max`. Same code path inside the helper.

- **Patch release `0.14.5`.** Fifth patch on the `0.14` line. CMake
  `VERSION` bumps from `0.14.4` to `0.14.5`; CHANGELOG gains a
  `0.14.5` block under "Fixed", "Performance", "Tests".

## Non-goals

- **No relocation of `requireSameGrid` to `core::detail::`.** Stays
  in `solve::detail::`; cross-namespace use from `fvm::` is the
  established pattern (`diff::detail::requireAxisExtent` is reused
  by `fvm::cellLaplacian` already).
- **No public-API change to `fvm::cellLaplacian` signature.** Same
  `cellLaplacian(psi, D)` interface; preconditions are added at the
  entry, not via new arguments.
- **No tightening beyond max-over-faces.** The "per-cell tight
  bound" alternative was considered and rejected — marginal extra
  perf gain over max-over-faces in most cases, requires changing
  the CFL helper signature from a scalar to a per-cell array.
- **No expansion of the FVM family.** No higher-order face-flux
  reconstruction, no rank-2 / vector FVM. Deferred items unchanged.
- **No new public symbols.** All changes are inside existing
  function bodies plus the (already-public) `Grid::operator==` from
  PR #20.

## Deferred questions

None. The audit thread starting from PR #19's stencil-aliasing
fix is complete:

- 0.14.1 (PR #19) — periodic stencil aliasing on small `N`.
- 0.14.2 (PR #20) — RHS-grid mismatch in `solve::integrate` /
  `axpyFresh`.
- 0.14.3 (PR #21) — odd-`N` off-by-one in spectral wavenumbers.
- 0.14.4 (PR #22) — audit of `spectral::partial` Nyquist-zeroing
  (no production-code change needed).
- **0.14.5 (this PR)** — D-grid mismatch + D > 0 contract + CFL
  tightening on the FVM heterogeneous-D path.
