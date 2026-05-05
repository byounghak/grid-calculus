# Plan: Audit `spectral::partial` / `spectral::gradient` Nyquist-zeroing path

## Group 1 — regression tests (the audit gate)

- New `test/spectral_partial_nyquist_test.cpp`. Each test constructs an input field with content at the highest positive harmonic on an odd-`N` axis, applies the spectral derivative, and compares to the closed-form analytical result.
  - **`spectral::partial<1, 0, 0>`** on `Nx ∈ {5, 7, 15, 31}` (rfft x-axis, first derivative). Input `sin(k_max · x)`, expected `cos(k_max · x) · k_max`.
  - **`spectral::partial<0, 1, 0>`** on `Ny ∈ {5, 7, 15, 31}` (full y-axis).
  - **`spectral::partial<0, 0, 1>`** on `Nz ∈ {5, 7, 15, 31}` (full z-axis).
  - **`spectral::partial<3, 0, 0>`**, **`<0, 3, 0>`**, **`<0, 0, 3>`** on the same odd-`N` axes (third derivative — odd per-axis count `3`). Expected output: `−cos(k_max · a) · k_max^3` for each axis `a`.
  - **`spectral::gradient`** on odd-`N` grids with input content at the highest positive harmonic on each axis (three sub-tests, one per active axis).
  - Tolerance: `< 1e-12` relative max-norm (spectral differentiation is exact for band-limited inputs to round-off; the multiplier residual at the highest harmonic is the gate).
- Wire into `test/CMakeLists.txt`.

**Commit message:** `test(spectral): regression tests for parity-gated Nyquist-zeroing on odd-N grids`

## Group 2 — Doxygen clarification

- Update the file-level Doxygen on [include/gridcalc/spectral/Derivatives.hpp](../../include/gridcalc/spectral/Derivatives.hpp) (lines 19-23) to make the parity gate explicit: "Nyquist mode is zeroed only when the axis extent is even AND the per-axis derivative count is odd. Odd `N` has no Nyquist mode by construction; the zeroing path is short-circuited and the highest positive harmonic at `(N − 1) / 2` is preserved with its full multiplier."
- Update the per-function Doxygen block on `spectral::partial` (the only function whose Doxygen explicitly mentions the Nyquist convention) for the same correction.

**Commit message:** `docs(spectral): clarify parity-gated Nyquist-zeroing in spectral derivatives`

## Group 3 — discharge the deferred follow-up

- Edit [specs/STATUS.md](../../specs/STATUS.md) — the "Decisions worth knowing" line 140 currently records "Nyquist mode is zeroed for odd-rank derivative axes" without naming the parity gate. Update to reflect what the audit confirmed: "the path is gated on `N % 2 == 0` AND odd per-axis derivative count; the new regression tests in [test/spectral_partial_nyquist_test.cpp](../../test/spectral_partial_nyquist_test.cpp) lock this." Add a "Phase 14 follow-up #4" entry under the post-Phase-14 follow-up list pointing to this spec dir + PR.
- Edit [specs/2026-05-04-fix-odd-N-fftfreq/requirements.md](../../specs/2026-05-04-fix-odd-N-fftfreq/requirements.md) "Deferred questions" section to mark the consumer-side audit resolved by this patch.

**Commit message:** `docs(specs): discharge spectral Nyquist-zeroing audit deferred from 0.14.3`

## Group 4 — version bump + CHANGELOG

- Bump `CMakeLists.txt` `VERSION` from `0.14.3` to `0.14.4`.
- Update `test/version_test.cpp` to assert `"0.14.4"`.
- Add a `## 0.14.4 — Audit: spectral::partial Nyquist-zeroing path on odd-N grids` block to `CHANGELOG.md` (above the `0.14.3` block) with `### Audited`, `### Tests`, `### Documentation` sections.

**Commit message:** `chore: release 0.14.4 -- audit-only patch closing the Nyquist-zeroing follow-up`
