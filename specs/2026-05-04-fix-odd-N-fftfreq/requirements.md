# Requirements: Fix odd-`N` off-by-one in `kyFull` / `kzFull`

## Bug report

`include/gridcalc/spectral/Wavenumbers.hpp:64,84` builds the
full-spectrum signed wavenumber index using

```cpp
const int n_signed = (n < Ny / 2) ? n : n - Ny;
```

The condition uses C++ integer division (truncating toward zero), so
for **odd `N`** the boundary lands one slot too low. For `Ny = 5`,
`Ny / 2 == 2`, and index `n = 2` (which `numpy.fft.fftfreq` reports
as `+2`, the highest positive harmonic) falls into the negative
branch and is mapped to `2 - 5 = -3` — a phantom mode at a frequency
that does not exist in the discrete spectrum.

Worked example with `Ny = 5`:

| `n` | `numpy.fftfreq` | current code | `2 * n < Ny` (correct) |
|---|---|---|---|
| 0 | `0`  | `0`  ✓ | `0`  ✓ |
| 1 | `+1` | `+1` ✓ | `+1` ✓ |
| 2 | `+2` | `−3` ✗ | `+2` ✓ |
| 3 | `−2` | `−2` ✓ | `−2` ✓ |
| 4 | `−1` | `−1` ✓ | `−1` ✓ |

**Impact.** Every spectral derivative on an odd-`Ny` or odd-`Nz`
grid uses the wrong wavenumber at index `(N − 1) / 2`. For
`spectral::laplacian` / `spectral::biharmonic` that mode carries the
largest `|k|^2` / `|k|^4` magnitude; the contaminated coefficient
injects a substantial wrong-frequency contribution into the result.

`kxRfft` is genuinely safe — the rfft half-spectrum holds only
non-negative harmonics, so the sign-flip path is absent. For odd
`Nx`, length `Nx/2 + 1 == (Nx − 1)/2 + 1`, and `kx[n] = scale * n`
correctly yields `[0, 1, ..., (Nx − 1)/2]`, matching
`numpy.fft.rfftfreq`.

## Why CI did not catch it

- `test/fd_fft_crosscheck_test.cpp:36` — sweep is
  `kSweepNs = {16, 32, 64, 128}`, all even.
- `test/spectral_basic_test.cpp` — uses `N = 32` and `N = 16`, both
  even.

The odd-`N` code path is uncovered. The spectral path is documented
as "verification only," but the FD–FFT cross-check would silently
disagree (or coincidentally pass at zero error if the test signal
has no content at index `(N − 1) / 2`) on odd-`N` grids.

## Decisions made for this feature

- **Fix: `(n < Ny / 2)` → `(2 * n < Ny)`** at lines 64 and 84.
  Equivalently `n < (Ny + 1) / 2`. The chosen form
  (`2 * n < Ny`) is the clearest "what numpy means by the boundary"
  reading and produces no integer-truncation footguns.
  - **Even-`N` regression check.** For even `Ny`, `Ny / 2 == (Ny + 1) / 2 - 1`
    and `2 * n < Ny` agrees with the original at every index
    (Nyquist `n = Ny/2` falls in the negative branch under both
    forms, mapping to `−Ny/2`, matching the standard
    `numpy.fft.fftfreq` convention). Existing even-`N` callers see
    byte-identical output.

- **Doxygen text correction.** The wrong formula `n_signed = n if n < N/2 else n - N`
  is also reproduced verbatim in the file-level block
  ([Wavenumbers.hpp:13](../../include/gridcalc/spectral/Wavenumbers.hpp#L13))
  and on `kyFull` / `kzFull`
  ([Wavenumbers.hpp:53](../../include/gridcalc/spectral/Wavenumbers.hpp#L53),
  [Wavenumbers.hpp:73](../../include/gridcalc/spectral/Wavenumbers.hpp#L73)).
  Update all three to state `2 * n < N` and call out the difference
  between even `N` (Nyquist mode at `n = N/2` lands in the negative
  branch by convention) and odd `N` (no Nyquist mode; highest
  positive harmonic at `n = (N − 1)/2`).

- **Consumer-side Nyquist-zeroing audit deferred to a follow-up.**
  `STATUS.md` line 140 records that `spectral::partial` "zeros the
  Nyquist mode for odd-rank derivative axes." If that path zeros by
  index unconditionally (e.g., `coeff[N/2] = 0`), it would silently
  zero a real positive-frequency mode at odd `N` (since `N/2` is the
  highest positive harmonic, not a Nyquist). This patch does **not**
  touch `spectral/Derivatives.hpp`; the audit is a separate
  follow-up. Recorded under "Deferred questions" so it is not
  dropped.

- **Test layout: three layers, defense-in-depth.**
  1. **Direct unit tests** in new `test/spectral_wavenumbers_test.cpp`:
     assert `kyFull(grid)` and `kzFull(grid)` match a hand-computed
     `numpy.fftfreq` reference array element-by-element at
     `N ∈ {4, 5, 6, 7, 15, 16, 31, 32}` (covering even / odd, small /
     large). Pin the exact bug.
  2. **FD–FFT cross-check sweep widened to include `N = 17`** in
     `test/fd_fft_crosscheck_test.cpp`. Adding one odd `N` to the
     existing `{16, 32, 64, 128}` sweep forces every existing slope-
     band assertion to pass on at least one odd-`N` grid. `N = 17`
     gives `h = 2π / 17 ≈ 0.37` — well inside the truncation-analysis
     regime, unlike 0.14.1's `N = 5` rejection (where `h ≈ 1.257`
     made the slope test unsound). The slope band stays
     `[Order − 0.5, Order + 0.5]`; the cross-check covers every
     existing FD↔FFT operator at every order.
  3. **Manufactured-solution highest-harmonic test** in the new
     wavenumbers test file: the strongest gate. Construct a field
     `psi(y) = sin(2π · ((Ny − 1)/2) · y / L)` for an odd `Ny`,
     apply `spectral::laplacian`, and assert the result equals the
     analytical eigenvalue `−|k|^2 · psi` to round-off. Without the
     fix this fails by orders of magnitude on odd-`N` grids
     (`|k|^2` at index `(Ny − 1)/2` becomes `|k|^2` at index
     `(Ny + 1)/2`, magnitudes differ); with the fix it passes to
     round-off. The same construction is exercised on z- and y-axis
     odd-`N` grids.

- **Patch release `0.14.3`.** Third patch release on the `0.14`
  line. CMake `VERSION` bumps from `0.14.2` to `0.14.3`; CHANGELOG
  gains a `0.14.3` block under "Fixed". `kyFull` / `kzFull`'s
  `\since` Doxygen tags note the fix
  (`\since 0.9.0 (function); 0.14.3 (odd-N off-by-one fix)`).

## Non-goals

- **No consumer-side audit** of `spectral/Derivatives.hpp`'s
  Nyquist-zeroing — split into a follow-up per the user's direction.
- **No `kxRfft` change.** It is correct for both parities by
  construction; left untouched (verified by reading the code:
  `kx[n] = scale * n` for `n ∈ [0, Nx/2 + 1)` matches
  `numpy.fft.rfftfreq` exactly).
- **No expansion of the FD–FFT cross-check sweep beyond one odd
  `N`.** Adding `N = 17` is the minimum that proves "odd-`N` is now
  exercised"; adding more (e.g., `N = 33`, `N = 65`) would multiply
  the fixture's runtime without adding coverage of any distinct
  failure mode.
- **No new public symbols.** The fix is a one-character change in
  each of two existing functions plus a Doxygen text correction; no
  new helper warranted, no new `\since` other than the in-text note.

## Deferred questions

- **Consumer-side Nyquist-zeroing audit (`spectral/Derivatives.hpp`).**
  Tracked as the next follow-up: read the operator that consumes
  `kyFull` / `kzFull`, verify whether the Nyquist-zero path is
  conditional on `N % 2 == 0`, and patch if not. Likely a separate
  `0.14.4` patch release.
- **User Guide / Developer Note "Nyquist mode" passages.** Phase 9's
  Developer Note chapter (Wavenumbers / Nyquist convention) may
  contain text that becomes inaccurate once the consumer-side audit
  lands. Tracked with the audit follow-up.
- **`Field` operator overloads, broader Field-equality work.** Out
  of scope; not motivated by this bug.
