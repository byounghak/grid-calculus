# Requirements: Fix silent stencil aliasing on small periodic axes

## Bug report

Centered finite-difference stencils of radius `r` are applied through the
input `Field`'s periodic index policy with no precondition that each axis
has at least `2r + 1` distinct samples. On any axis with `N ≤ 2r`, the
periodic wrap aliases two or more offsets in `[-r, r]` to the same grid
point, and the advertised stencil silently degrades to a different (and
generally non-zero, non-truncation-controlled) operator.

Concrete failure case from the bug report: `laplacian<4>(field)` on a
grid with `Nx = 4`. Offsets `[-2, -1, 0, 1, 2]` wrap to indices
`{2, 3, 0, 1, 2}` — index `2` appears twice. The two `-1/12` weights at
offsets ±2 sum to `-1/6` at the same grid sample, so the `O(h^4)`
five-point stencil collapses to a four-point stencil with weights
`{4/3, -8/3, 4/3, -1/6}` (after summation) — the published convergence
order no longer holds, and there is no compile-time, link-time, or
runtime indication that the user has stepped outside the operator's
domain of validity.

## Affected operators

| Stencil family                 | Order | radius | Min `N` (per active axis) |
|--------------------------------|-------|--------|---------------------------|
| `Coefficients<2>` (Laplacian)  | 2     | 1      | 3                         |
| `FirstDerivative<2>`           | 2     | 1      | 3                         |
| `Coefficients<4>` (Laplacian)  | 4     | 2      | 5                         |
| `FirstDerivative<4>`           | 4     | 2      | 5                         |
| `ThirdDerivative<2>`           | 2     | 2      | 5                         |
| `FourthDerivative<2>`          | 2     | 2      | 5                         |
| `ThirdDerivative<4>`           | 4     | 3      | 7                         |
| `FourthDerivative<4>`          | 4     | 3      | 7                         |

The bug touches every Phase 7 and Phase 8 high-order operator
(`laplacian<4>`, `gradient<4>`, `divergence<4>`, all 31 d-prefix mixed
partials at both orders, `biharmonic<2>` and `biharmonic<4>`) plus the
Phase 13 rank-raising `gradient` / `divergence` overloads that inherit
the same `Order` parameter. `Coefficients<2>` and `FirstDerivative<2>`
fail only on `N ∈ {1, 2}` — essentially degenerate, but the precondition
is added uniformly for consistency.

`fvm::cellLaplacian` (Phase 14) uses radius-1 face fluxes and is
therefore safe for any `N ≥ 3` the `Grid` constructor admits. The
precondition is added defensively to keep the operator-entry contract
uniform across `diff::` and `fvm::`.

## Why CI did not catch it

The Phase 9 FD–FFT cross-check fixture (`test/fd_fft_crosscheck_test.cpp`)
parameterizes over the grid sweep `N ∈ {16, 32, 64, 128}` (recorded in
`specs/STATUS.md` line 139). All four values exceed every operator's
`2r + 1` threshold, so the alias never triggers in CI. No prior test
exercises the smallest valid `N` for a radius-2 or radius-3 stencil.

## Decisions made for this feature

- **Failure mode: `throw std::invalid_argument`** at the operator entry
  point if any axis the operator reads has `2 <= N <= 2 * radius`.
  Matches the Grid-constructor pattern
  (`include/gridcalc/core/Grid.hpp`); fires identically in Debug and
  Release; the production / industrial mission target makes "silent in
  Release" unacceptable. The error message names the offending axis,
  its `N`, the stencil's radius, and the required minimum, so the user
  can fix the call site without reading source. One branch per
  operator call (not per grid point), so the happy-path cost is
  negligible.

- **Degenerate-axis carve-out: `N = 1` is accepted, regardless of
  stencil radius.** Mathematical justification: at `N = 1` every
  stencil offset `[-r, r]` wraps to index `0`, all `2 * radius + 1`
  reads return the single cell value `psi(0)`, and the centered FD
  weight tables in this codebase (1st through 4th derivative) sum to
  zero by construction (the n-th derivative of a constant is zero).
  The aliased result is therefore *exactly* `0`, which is the correct
  value for a degenerate single-cell periodic axis where there is no
  spatial variation. Practically, this preserves the existing 1D-style
  test pattern `Grid g(4, 1, 1, ...)` used by Phase 14's FVM unit
  tests. Only the *partial-alias* range `[2, 2 * radius]` is
  rejected — that is the genuine "silent-degradation" zone the bug
  report is concerned with. The carve-out is documented inline in the
  helper's Doxygen and exercised by `DegenerateAxis_*_AcceptsN1` tests
  in `test/stencil_axis_extent_test.cpp`.
- **Single helper.** The check is a free function
  `gridcalc::diff::detail::requireAxisExtent(axis_label, N, radius)` in
  a new header `include/gridcalc/diff/detail/PreconditionAxisExtent.hpp`,
  reused by every operator. Keeps the error-message format consistent
  and the per-call-site footprint to one line.
- **Per-axis active set.** For Phase 8 mixed partials accessed through
  `diff::detail::mixedDerivative<Nx, Ny, Nz, Order>`, the check fires
  only on axes with the corresponding per-axis derivative count `> 0`
  (an axis with `n_a = 0` contributes a single delta-stencil at offset
  `0`, which never reads off-axis samples). For all other operators
  (Laplacian, gradient, divergence, biharmonic), every axis is active.
- **Test scope: throw per radius tier + smallest-valid-N boundary
  (no-throw + finiteness gate).** Throw tests at `N = 2 * radius` (the
  largest aliasing case) for each radius tier (`1, 2, 3`) across every
  affected operator. In addition, smallest-valid-N tests at
  `N = 2 * radius + 1` for every Phase 7 / Phase 8 / Phase 13 /
  Phase 14 operator, asserting no-throw and finite output (no NaN /
  Inf). Both groups live in the new `test/stencil_axis_extent_test.cpp`.

  **Why a no-throw + finiteness gate, not a slope-sweep extension.**
  Extending the Phase 9 FD–FFT cross-check fixture's grid sweep
  (`{16, 32, 64, 128}`) to include `N = 2 * radius + 1` would be unsound:
  at `N = 5` the grid spacing `h = 2π / 5 ≈ 1.257` is large enough that
  the truncation-error analysis underpinning the slope test breaks down,
  and the slope-band assertion `[Order - 0.5, Order + 0.5]` would be
  polluted by the small-N point. The boundary check captures the actual
  failure mode under the bug (silent aliasing produces non-finite or
  catastrophically-wrong outputs) without claiming a truncation order
  the test cannot honestly enforce. Numerical correctness at
  production-relevant `N` is already covered by the existing FD–FFT
  sweep at the unchanged `{16, 32, 64, 128}` entries.
- **`fvm::cellLaplacian` gets the same check.** Defense-in-depth; the
  FVM cell-flux discretization at Phase 14 ships with radius `1`, so
  the check is currently a no-op for any `N ≥ 3`. Keeps the entry
  contract uniform across `diff::` and `fvm::`; if a future
  higher-order FVM face-flux variant lands (radius `≥ 2`), the check is
  already in place at the right call site.
- **Patch release `0.14.1`.** First patch release for this project.
  CMake `VERSION` bumps from `0.14.0` to `0.14.1`; CHANGELOG gains a
  `0.14.1` block under "Fixed" describing the silent-aliasing bug and
  the precondition. The `\since` Doxygen tags on the new helper carry
  `\since 0.14.1`; existing operator Doxygen blocks gain a `\throws`
  line referencing the precondition.

## Non-goals

- **No spectral check.** `gridcalc::spectral::*` operates on global
  Fourier coefficients with no per-axis stencil radius. The aliasing
  bug does not apply; spectral operators are unmodified.
- **No public type for the helper.** `requireAxisExtent` lives under
  `gridcalc::diff::detail` (and is re-used by `fvm::cellLaplacian` via
  a relative include — `fvm` is the only external consumer for now).
  The helper is an implementation detail of the precondition; it is not
  documented for end users. The Phase 10 `\since`-tag lint regex covers
  the public-symbol set; the helper carries Doxygen for engineers
  reading the codebase but is excluded from the public API surface.
- **No global "max radius" lookup on `Grid`.** A `Grid::maxRadius()`
  method was considered as a way to hoist the check upstream of the
  operator entry. Rejected: the radius is a property of the *operator*,
  not the *grid*, and storing it on `Grid` would couple the storage
  type to the differentiation surface in a way that breaks the
  Phase 1 / Phase 2 separation (cf. `core/Grid.hpp` docstring).
- **No widening of the FD–FFT cross-check N sweep.** Adding small-N
  entries (`N = 5, 7, ...`) was considered and rejected — see the
  no-throw + finiteness rationale under "Decisions made for this
  feature." The existing `{16, 32, 64, 128}` sweep is unchanged.
- **No CHANGELOG migration of historical entries.** The `0.14.0` block
  stays as-is; this is a clean additive `0.14.1` block.

## Deferred questions

- **Higher-order FVM.** If a Phase 21+ FVM variant introduces a
  radius-2 face-flux reconstruction (MUSCL etc.), the precondition is
  already wired through `fvm::cellLaplacian`. Tracked as a deferred
  item in `STATUS.md`'s "Open / deferred items" section already.
- **Compile-time radius check.** A `static_assert` would be ideal but
  `Grid` extents are runtime values from `int` constructor arguments;
  no compile-time path is currently feasible. If a future
  `FixedExtentGrid<Nx, Ny, Nz>` lands (none planned), it could carry a
  compile-time form of the same check.
