# Validation: Fix odd-`N` off-by-one in `kyFull` / `kzFull`

## Build

- [ ] CI green on Ubuntu GCC.
- [ ] CI green on Ubuntu Clang.
- [ ] CI green on the `clang-debug-nofft` configuration. (This patch is gated on `GRIDCALC_USE_FFT`; the no-FFT build skips the changed code path entirely. The build must still configure / link / test cleanly.)
- [ ] Warnings clean under the project's strict warning set.
- [ ] `clang-format` and `clang-tidy` pass.
- [ ] `docs-build` job green.
- [ ] `docs-lint` job green. `kyFull` / `kzFull`'s `\since` tags note the precondition (`\since 0.9.0 (function); 0.14.3 (odd-N off-by-one fix)`).

## Tests

- [ ] `test/spectral_wavenumbers_test.cpp` exists and contains:
  - [ ] Direct array-equality tests for `kyFull` at `Ny ∈ {4, 5, 6, 7, 15, 16, 31, 32}` against a hand-computed numpy `fftfreq` reference.
  - [ ] Direct array-equality tests for `kzFull` at the same `Nz` values.
  - [ ] `kxRfft` regression at odd `Nx ∈ {5, 7}` (no behavior change expected; documents the contract).
  - [ ] Manufactured-solution highest-positive-harmonic test on a y-axis odd-`Ny` grid: `psi(y) = sin(k_max · y)` recovers `−k_max² · psi` from `spectral::laplacian` to round-off.
  - [ ] Same on a z-axis odd-`Nz` grid.
- [ ] `test/fd_fft_crosscheck_test.cpp` — `kSweepNs` widened to include `N = 17`; existing `{16, 32, 64, 128}` entries unchanged. Every existing FD↔FFT operator slope-band assertion (`[Order − 0.5, Order + 0.5]`) passes on the widened sweep.
- [ ] Full test suite passes on `clang-debug` (existing 316 tests + new tests).
- [ ] Full test suite passes on `clang-debug-nofft` (existing 240 tests; the spectral path is `#ifdef`-skipped, so the new `spectral_wavenumbers_test.cpp` adds zero tests under this config).

## Documentation

- [ ] Wavenumbers.hpp file-level Doxygen block (line 13) updated: the wrong formula `n_signed = n if n < N/2 else n - N` is replaced with `2 * n < N` and the even/odd boundary distinction is called out.
- [ ] Per-function Doxygen on `kyFull` (line 53) and `kzFull` (line 73) is updated with the same correction.
- [ ] `\since` tags note the precondition addition.
- [ ] `CHANGELOG.md` carries a new `0.14.3` block above the existing `0.14.2` block with `### Fixed`, `### Tests`, `### Documentation` sections.
- [ ] `specs/STATUS.md` "Last updated" line bumped; "Phase 14 follow-up #3" entry added.

## Versioning

- [ ] `CMakeLists.txt` project `VERSION` is `0.14.3`.
- [ ] `test/version_test.cpp` asserts `"0.14.3"`.

## Manual verification

- [ ] On a fresh build, `kyFull(Grid(8, 5, 8, ...))` produces `[0, 1, 2, -2, -1] * (2π / (5 * hy))` element-wise (the numpy `fftfreq` reference for `Ny = 5`).
- [ ] On the same fresh build, `spectral::laplacian` on `sin(k_max · y)` at `Ny = 5` recovers `−k_max² · sin(k_max · y)` to round-off; before the fix this fails by O(|k_max|²) magnitude.
- [ ] The existing `clang-debug` test surface (316 tests as of 0.14.2) is unaffected by the fix on any even-`N` grid (regression check via the unchanged `{16, 32, 64, 128}` slope sweep).

## Deferred for follow-up

- [ ] Audit `spectral/Derivatives.hpp` (`spectral::partial`'s "Nyquist mode is zeroed for odd-rank derivative axes" path per `STATUS.md` line 140) — verify the zeroing is conditional on `N % 2 == 0`, patch if not. Tracked as the next follow-up; not in this PR.
