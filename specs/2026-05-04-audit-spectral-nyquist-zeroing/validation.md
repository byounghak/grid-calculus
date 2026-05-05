# Validation: Audit `spectral::partial` / `spectral::gradient` Nyquist-zeroing path

## Build

- [ ] CI green on Ubuntu GCC.
- [ ] CI green on Ubuntu Clang.
- [ ] CI green on the `clang-debug-nofft` configuration. (This audit's tests are gated on `GRIDCALC_USE_FFT`; the no-FFT build skips the changed code path. The build must still configure / link / test cleanly.)
- [ ] Warnings clean under the project's strict warning set.
- [ ] `clang-format` and `clang-tidy` pass.
- [ ] `docs-build` job green.
- [ ] `docs-lint` job green.

## Tests

- [ ] `test/spectral_partial_nyquist_test.cpp` exists and exercises:
  - [ ] `spectral::partial<1, 0, 0>` on `Nx ∈ {5, 7, 15, 31}` with content at the highest positive harmonic; result matches `cos(k_max · x) · k_max` to round-off.
  - [ ] `spectral::partial<0, 1, 0>` on `Ny ∈ {5, 7, 15, 31}` analogously.
  - [ ] `spectral::partial<0, 0, 1>` on `Nz ∈ {5, 7, 15, 31}` analogously.
  - [ ] `spectral::partial<3, 0, 0>`, `<0, 3, 0>`, `<0, 0, 3>` on the same odd-`N` axes; result matches `−cos(k_max · a) · k_max^3` to round-off.
  - [ ] `spectral::gradient` on odd-`N` grids with content at the highest positive harmonic on each axis (three sub-tests, one per active axis); the relevant component matches the analytical first-derivative.
- [ ] Full test suite passes on `clang-debug` (existing 347 tests + new tests).
- [ ] Full test suite passes on `clang-debug-nofft` (existing 240 tests; the new test file is `GRIDCALC_USE_FFT`-gated).

## Documentation

- [ ] [include/gridcalc/spectral/Derivatives.hpp](../../include/gridcalc/spectral/Derivatives.hpp) file-level Doxygen block (lines 19-23) updated to make the parity gate explicit.
- [ ] Per-function Doxygen on `spectral::partial` updated to note the parity gate.
- [ ] No `\since` tag changes (no new public symbols).
- [ ] `CHANGELOG.md` carries a new `0.14.4` block above the existing `0.14.3` block with `### Audited`, `### Tests`, `### Documentation` sections.
- [ ] `specs/STATUS.md` "Decisions worth knowing" line 140 updated to reflect the audit conclusion; "Phase 14 follow-up #4" entry added.
- [ ] `specs/2026-05-04-fix-odd-N-fftfreq/requirements.md` "Deferred questions" section updated to mark the consumer-side audit resolved.

## Versioning

- [ ] `CMakeLists.txt` project `VERSION` is `0.14.4`.
- [ ] `test/version_test.cpp` asserts `"0.14.4"`.

## Manual verification

- [ ] On a fresh build, `spectral::partial<0, 1, 0>(sin(k_max · y))` on `Ny = 5` recovers `cos(k_max · y) · k_max` to round-off (`< 1e-12` relative). This is the "would-fail-without-correct-parity-gating" gate that pins the audit.
- [ ] No production-code changes to `Derivatives.hpp`'s zeroing logic (audit-only patch).
