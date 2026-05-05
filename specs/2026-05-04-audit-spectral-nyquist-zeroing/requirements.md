# Requirements: Audit `spectral::partial` / `spectral::gradient` Nyquist-zeroing path

## Audit motivation

The 0.14.3 patch (`fix: off-by-one in spectral wavenumbers on odd-N grids`)
recorded a deferred follow-up: "Audit `spectral::partial`'s 'Nyquist
mode is zeroed for odd-rank derivative axes' path — if zeroed by
index unconditionally rather than gated on `N % 2 == 0`, it would
silently zero a real positive-frequency mode at odd `N` (since
`N / 2` is the highest positive harmonic, not a Nyquist)."
Recorded under `STATUS.md` "Decisions worth knowing" and in the
0.14.3 spec dir's `Deferred questions` section.

This patch discharges that audit.

## Audit conclusion

**The consumer code is already correct.** Reading
[include/gridcalc/spectral/Derivatives.hpp](../../include/gridcalc/spectral/Derivatives.hpp):

- [`spectral::partial<Nx, Ny, Nz>`](../../include/gridcalc/spectral/Derivatives.hpp#L103-L117)
  zeros the Nyquist mode along axis `a` only when both
  `Nx_grid % 2 == 0` (axis is even) **and** the axis derivative
  count is odd (`if constexpr (Nx % 2 == 1)`). For odd `N` the
  parity check short-circuits and no zeroing occurs.
- [`spectral::gradient`](../../include/gridcalc/spectral/Derivatives.hpp#L221-L229)
  applies the same parity-gated check on each Cartesian component.
- [`spectral::laplacian`](../../include/gridcalc/spectral/Derivatives.hpp#L135)
  and [`spectral::biharmonic`](../../include/gridcalc/spectral/Derivatives.hpp#L166)
  do no Nyquist zeroing at all — their per-axis derivative counts
  are even (`(2, 0, 0) + (0, 2, 0) + (0, 0, 2)` for laplacian and
  the 4th-order analogue for biharmonic), so the `(±π/h)^N` sign
  ambiguity does not arise.

The Phase 9 FD–FFT cross-check sweep at `N = 17` (added in 0.14.3)
exercises every FD↔FFT operator at every order on at least one
odd-`N` grid and passes — corroborating that the consumer code is
correct under odd `N` for the operators that have content at
non-degenerate harmonics.

What the existing test surface does **not** cover: an input with
content *at the highest positive harmonic* on an odd-`N` axis,
combined with an *odd* per-axis derivative count. The 0.14.3
manufactured-solution highest-harmonic test exercises this for
`spectral::laplacian` (even derivative count, no zeroing path),
not for `spectral::partial<0, 1, 0>` / `<3, 0, 0>` / etc. The new
regression tests in this patch close that gap.

## Decisions made for this feature

- **Verify-and-discharge approach.** No production-code change to
  `Derivatives.hpp`'s zeroing logic — it is correct. Patch ships:
  (a) regression tests that pin the correctness, (b) a Doxygen
  clarification on the file-level docstring making the parity-gated
  rule explicit, and (c) a discharge note in `STATUS.md` and the
  0.14.3 spec dir's deferred-questions section pointing to this
  audit.

- **Test scope.** New `test/spectral_partial_nyquist_test.cpp`
  exercises every "odd-N grid + odd per-axis derivative count +
  input content at the highest positive harmonic" combination that
  could plausibly hide a bug:
  - **`spectral::partial<1, 0, 0>`** (rfft x-axis, first derivative)
    on `Nx ∈ {5, 7, 15, 31}` — input `sin(k_max · x)`,
    `k_max = 2π · ((Nx − 1)/2) / L`; expected output `cos(k_max · x) · k_max`.
  - **`spectral::partial<0, 1, 0>`** (full y-axis) and
    **`spectral::partial<0, 0, 1>`** (full z-axis) on
    `Ny / Nz ∈ {5, 7, 15, 31}`.
  - **`spectral::partial<3, 0, 0>`**, **`<0, 3, 0>`**, **`<0, 0, 3>`**
    on the same odd-`N` axes (third derivative — odd count `3`).
  - **`spectral::gradient`** on odd-`N` grids with input content at
    the highest positive harmonic on each axis; the relevant
    component must recover the analytical result.

- **No production-code change to the zeroing logic.** Confirmed
  correct by reading + by the new regression tests passing without
  any change.

- **Doxygen clarification.** [Derivatives.hpp:19-23](../../include/gridcalc/spectral/Derivatives.hpp#L19)
  currently states "for any axis with an odd per-axis derivative
  count, the mode at `n_a = N_a / 2` (Nyquist) is set to zero"
  without naming the parity gate. Update to: "for any axis with an
  odd per-axis derivative count *and an even axis extent `N_a`*,
  the mode at `n_a = N_a / 2` (the Nyquist mode for even `N_a`) is
  set to zero. Odd `N_a` has no Nyquist mode by construction; the
  zeroing path is short-circuited and the highest positive harmonic
  at `(N_a − 1) / 2` is preserved with its full multiplier." Same
  clarification in the per-function `partial` Doxygen block.

- **Discharge mechanics.** Update both `STATUS.md` "Decisions worth
  knowing" line 140 and the 0.14.3 spec dir's
  `Deferred questions` section to mark the audit complete and
  point readers to this spec dir + the new regression test file
  as the gate.

- **Patch release `0.14.4`.** Audit-only patch but ships test code,
  Doxygen content, and a CHANGELOG entry; the version bump is
  warranted by the project's per-PR version-bump pattern. CMake
  `VERSION` bumps from `0.14.3` to `0.14.4`; CHANGELOG block under
  "Audited" / "Tests" / "Documentation" sections.

## Non-goals

- **No code change to `Derivatives.hpp`'s zeroing logic.** Already
  correct.
- **No FD–FFT cross-check sweep changes.** The 0.14.3 widening to
  `{16, 17, 32, 64, 128}` is sufficient; this audit's gate is the
  new dedicated regression test file.
- **No revisit of the 0.14.3 fix to `kyFull` / `kzFull`.** Already
  landed; this audit is on the consumer side.
- **No broader audit of `Fft.hpp` / `spectrumIndex` / other spectral
  helpers.** Out of scope; the deferred-follow-up was specifically
  about the Nyquist-zeroing path.

## Deferred questions

None remaining from the 0.14.3 follow-up. If a future audit surfaces
something in `Fft.hpp` or `spectrumIndex`, that's a separate item.
