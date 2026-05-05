# Plan: Fix odd-`N` off-by-one in `kyFull` / `kzFull`

## Group 1 — patch the off-by-one + correct Doxygen

- Edit [include/gridcalc/spectral/Wavenumbers.hpp](../../include/gridcalc/spectral/Wavenumbers.hpp):
  - Line 64 (`kyFull`): `(n < Ny / 2)` → `(2 * n < Ny)`.
  - Line 84 (`kzFull`): `(n < Nz / 2)` → `(2 * n < Nz)`.
  - File-level Doxygen (line 13): correct the formula and call out the even/odd boundary distinction.
  - Per-function Doxygen on `kyFull` (line 53) and `kzFull` (line 73): same correction.
  - `\since` tags on `kyFull` and `kzFull`: append `; 0.14.3 (odd-N off-by-one fix)`.

**Commit message:** `fix(spectral): off-by-one in kyFull / kzFull on odd-N grids`

## Group 2 — direct unit tests + widen FD–FFT cross-check sweep

- New `test/spectral_wavenumbers_test.cpp`:
  - **Direct array-equality** sub-fixture: assert `kyFull(grid)` and `kzFull(grid)` match a hand-computed numpy `fftfreq` reference array element-by-element. One test per parity tier:
    - Even: `Ny ∈ {4, 6, 16, 32}`, `Nz` symmetrically.
    - Odd: `Ny ∈ {5, 7, 15, 31}`, `Nz` symmetrically.
  - **`kxRfft` regression** — confirm the rfft path is unchanged at odd `Nx` (nothing to fix; document the contract via test).
- Edit [test/fd_fft_crosscheck_test.cpp](../../test/fd_fft_crosscheck_test.cpp): widen `kSweepNs` from `{16, 32, 64, 128}` to `{16, 17, 32, 64, 128}`. Adding `N = 17` is the minimum increment that exercises the odd-`N` code path; `h = 2π/17 ≈ 0.37` keeps the truncation-analysis regime intact (unlike 0.14.1's rejected `N = 5` at `h ≈ 1.257`). Slope band assertion stays `[Order − 0.5, Order + 0.5]`; every existing FD↔FFT operator at every order is exercised on at least one odd-`N` grid.

**Commit message:** `test(spectral): cover odd-N wavenumbers + add N=17 to FD-FFT sweep`

## Group 3 — manufactured-solution highest-harmonic acceptance test

- Append to `test/spectral_wavenumbers_test.cpp`: a manufactured-solution acceptance test exercising the highest-positive-harmonic mode at odd `N` — the exact mode that flipped under the bug. For an odd `Ny`, `psi(x, y, z) = sin(k_max · y)` with `k_max = 2π * (Ny − 1)/2 / L`; `spectral::laplacian(psi)` must equal `−k_max^2 · psi` to round-off. Repeat for the y- and z-axis on odd grids. Without the fix the test fails by `O(|k_max|^2)` magnitude (the mode flipped to a different frequency); with the fix it passes to round-off.

**Commit message:** `test(spectral): add highest-positive-harmonic acceptance test on odd N`

## Group 4 — version bump + CHANGELOG + STATUS

- Bump `CMakeLists.txt` `VERSION` from `0.14.2` to `0.14.3`.
- Update `test/version_test.cpp` to assert `"0.14.3"`.
- Add a `## 0.14.3 — Fix: off-by-one in spectral wavenumbers on odd-N grids` block to `CHANGELOG.md` (above the `0.14.2` block) with `### Fixed`, `### Tests`, `### Documentation` sections.
- Update `specs/STATUS.md`: add a "Phase 14 follow-up #3" entry with the merge-commit hash to be filled in post-merge; bump the **Last updated** line.

**Commit message:** `chore: release 0.14.3 -- odd-N fftfreq off-by-one fix`
