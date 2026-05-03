# Plan: Phase 9 — Spectral verification harness

## Group 1 — Vendor PocketFFT and wire into the build

- Create `extern/pocketfft/` and commit `pocketfft_hdronly.h` from upstream `mreineck/pocketfft` at a specific revision. Add `extern/pocketfft/VERSION` recording the upstream URL, revision SHA, and the date the vendor snapshot was taken.
- Add a CMake option `option(GRIDCALC_USE_FFT "Build the spectral verification harness (PocketFFT-backed)" ON)`.
- When `GRIDCALC_USE_FFT` is on:
  - `target_include_directories(gridcalc INTERFACE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/extern/pocketfft> ...)`.
  - `target_compile_definitions(gridcalc INTERFACE GRIDCALC_USE_FFT)`.
- Update `specs/tech-stack.md`'s Documentation / Dependencies section to record the vendored revision and the build flag.

**Commit message:** `build: vendor PocketFFT header and add GRIDCALC_USE_FFT build flag`

## Group 2 — Wavenumber helper

- Create `include/gridcalc/spectral/Wavenumbers.hpp`:
  - `core::Field<double>`-shaped wavenumber arrays per axis. Specifically, three `std::array`-like helpers that compute $k_a[n_a]$:
    - X (rfft half-spectrum): `kxRfft(grid)` returns a `std::vector<double>` of length `Nx/2 + 1` with $k_x[n_x] = 2\pi n_x / (N_x h_x)$ for $n_x \in [0, N_x/2]$.
    - Y / Z (full spectra): `kyFull(grid)`, `kzFull(grid)` of length `Ny`, `Nz`. Uses standard `numpy.fft`-style sign convention: $k_a[n_a] = 2\pi n_a^{\text{signed}} / (N_a h_a)$ where $n_a^{\text{signed}} = n_a$ if $n_a < N_a/2$ else $n_a - N_a$.
  - Doxygen blocks with `\since 0.9.0`.
- All three helpers are non-member functions taking `const core::Grid&`. They allocate small per-axis vectors and return them by value (cheap; the operator headers cache them once per call).

**Commit message:** `feat(spectral): add Wavenumbers.hpp helper for FFT k-space conventions`

## Group 3 — `spectral::Fft` thin wrapper

- Create `include/gridcalc/spectral/Fft.hpp`:
  - `spectral::ComplexField` alias = `std::vector<std::complex<double>>` of size `(Nx/2+1) * Ny * Nz`, sharing the i-fastest layout convention (see `STATUS.md` "i-fastest" note).
  - `spectral::forwardR2C(const core::Field<double>& field) -> ComplexField` — forward 3D rfft via `pocketfft::r2c`. **No** normalization (the standard is to divide by $N_x N_y N_z$ on the inverse transform).
  - `spectral::backwardC2R(const ComplexField& spectrum, const core::Grid& grid) -> core::Field<double>` — backward 3D crfft via `pocketfft::c2r`. Divides the result by $N_x N_y N_z$.
  - Header includes `<extern/pocketfft/pocketfft_hdronly.h>` (relative include via the configured `target_include_directories`).
  - `\since 0.9.0`.
- Whole header is gated by `#ifdef GRIDCALC_USE_FFT`. If the flag is off, the header is empty (or absent from `target_include_directories`).

**Commit message:** `feat(spectral): add Fft.hpp r2c / c2r wrapper over PocketFFT`

## Group 4 — Spectral derivative operators

- Create `include/gridcalc/spectral/Derivatives.hpp`:
  - `spectral::partial<int Nx, int Ny, int Nz>(const core::Field<double>& field) -> core::Field<double>`:
    1. Forward rfft → `ComplexField`.
    2. Build per-axis $k_a$ arrays via `Wavenumbers.hpp`.
    3. For each spectrum element $(n_x, n_y, n_z)$, multiply by $(i k_x)^{N_x} (i k_y)^{N_y} (i k_z)^{N_z}$.
    4. **Nyquist zeroing**: for each axis $a$ with $N_a$ odd, zero the spectrum at the Nyquist index ($n_a = N_a/2$ on the relevant axis). For x: $n_x = N_x/2$ (last index of half-spectrum). For y/z: $n_a = N_a/2$ if $N_a$ even.
    5. Backward c2r → real `Field<double>`.
  - Named convenience wrappers, all 1–2 line forwarders:
    - `spectral::laplacian(field)` = `spectral::partial<2,0,0>(field) + partial<0,2,0>(field) + partial<0,0,2>(field)` — but implemented in **one fused pass** that multiplies the spectrum by $-|\mathbf{k}|^2$ once (single forward/backward FFT pair).
    - `spectral::biharmonic(field)` — same, multiplies by $|\mathbf{k}|^4$ in one pass.
    - `spectral::gradient(field) -> core::Field<core::Vec3d>` — three spectral first-derivatives. Computed in a single forward FFT and three backward FFTs (one per Cartesian component); results are packed into the `Vec3d` field.
  - Doxygen blocks with `\since 0.9.0`.
- All declarations / definitions gated by `#ifdef GRIDCALC_USE_FFT`.

**Commit message:** `feat(spectral): add partial<Nx,Ny,Nz>, laplacian, biharmonic, gradient`

## Group 5 — FD–FFT cross-check fixture

- Add `test/fd_fft_crosscheck_test.cpp`:
  - Helper template `crossCheckSlope<int Nx, int Ny, int Nz, int Order>(applierFD)`:
    - For each `N ∈ {16, 32, 64, 128}`:
      - Build `Grid(N, N, N, h, h, h)` with `h = 2π/N`.
      - Build `Field<double>` $\psi$ from $\sin(x + y + z)$ — `k = (1, 1, 1)`, well below Nyquist on every grid.
      - Apply the FD operator (`applierFD(psi)`) → `Field<double>`.
      - Apply `spectral::partial<Nx, Ny, Nz>(psi)` → `Field<double>`.
      - Compute `maxAbsDiff / max(maxAbs(spectral), 1e-300)`.
    - Return log-log slope of the four discrepancies vs $h$.
  - Macro `EXPECT_FD_FFT_SLOPE(NAME_CHAIN, NX, NY, NZ, ORDER)` that:
    - Generates `TEST(FdFftCrossCheck, NAME_CHAIN ## _Order ## ORDER)`.
    - Calls `crossCheckSlope<NX, NY, NZ, ORDER>([](const Field<double>& f){ return NAME_CHAIN<ORDER>(f); })` and asserts the slope lies in `[Order - 0.5, Order + 0.5]`.
  - Phase-1 / -2 special cases that don't fit the `partial<Nx, Ny, Nz>` template directly (gradient returns `Field<Vec3d>`, divergence consumes `Field<Vec3d>`):
    - For `gradient`: cross-check each component against `spectral::partial<1,0,0>`, `<0,1,0>`, `<0,0,1>` separately.
    - For `divergence`: build a `Field<Vec3d>` source and cross-check `divergence<Order>` against `spectral::partial<1,0,0> + partial<0,1,0> + partial<0,0,1>` applied component-wise to the `Vec3d` field.
  - For `biharmonic` (Phase 8) — cross-check vs `spectral::biharmonic`.
  - 70 macro invocations covering Phases 1, 2, 8.
- Wire `test/fd_fft_crosscheck_test.cpp` into `test/CMakeLists.txt` under an `if(GRIDCALC_USE_FFT)` guard.
- Add a small `test/spectral_basic_test.cpp` smoke-testing the FFT round-trip and the named wrappers `spectral::{laplacian, biharmonic, gradient}` against $\sin(\mathbf{k}\cdot\mathbf{x})$ at single $N = 32$ (sanity gate before the heavy cross-check fixture).

**Commit message:** `test(spectral): FD-FFT cross-check fixture covering all Phase 1-8 operators`

## Group 6 — Phase 9 doc-notes

- `docs/user-guide/notes/phase-9-spectral-verification.md`:
  - The hybrid `spectral::partial<Nx, Ny, Nz>(field)` + named wrappers, when to use each.
  - **Strong "verification only" admonition**: per `mission.md` and the roadmap, FD operators are the production code; `spectral::*` exists to gate FD changes. Do not use `spectral::*` in production solvers.
  - Build flag `-DGRIDCALC_USE_FFT={ON,OFF}`.
  - Worked example: a CI-style assertion harness that runs `spectral::laplacian` and `diff::laplacian<4>` on the same input, computes the discrepancy, and asserts an order-4 bound.
- `docs/developer-note/notes/phase-9-spectral-verification.md` — five-section structure per `specs/CLAUDE.md` step 4d:
  1. **Theory.** The spectral derivative on a periodic grid as multiplication by $i\mathbf{k}$ in Fourier space; equivalence to the discrete derivative on band-limited inputs; why the spectral path is **not** equivalent to the FD path on non-band-limited inputs (aliasing). Wavenumber sign convention. Nyquist-mode aliasing for odd-rank derivatives — formal statement and remedy.
  2. **Math derivation.** Why $(\partial^{\boldsymbol\alpha}\psi)\widehat{}(\mathbf{k}) = (i\mathbf{k})^{\boldsymbol\alpha}\widehat\psi(\mathbf{k})$ for periodic $\psi$. Why FD–FFT discrepancy on a single trig mode below Nyquist scales as $\mathcal{O}(h^p)$ where $p$ is the FD order.
  3. **Algorithm.** PocketFFT's r2c / c2r interface; the half-spectrum layout; per-axis wavenumber construction; the `spectral::partial<Nx, Ny, Nz>` template; named wrappers as fused single-FFT-pair specializations; cross-check fixture structure. File-path references.
  4. **Design decisions.** Vendored PocketFFT vs `FetchContent`; hybrid API vs full named family; rfft vs cfft; $N \in \{16, 32, 64, 128\}$ sweep grids; verification-only positioning. Cite `requirements.md`.
  5. **References.** PocketFFT upstream; Trefethen *Spectral Methods in MATLAB* (Chapter 3, FFT-based differentiation, ISBN 978-0-89871-465-4 + DOI); Frigo–Johnson *FFTW* paper for the canonical real-input wavenumber-layout convention (DOI); Boyd *Chebyshev and Fourier Spectral Methods* for the band-limited-input requirement.

**Commit message:** `docs(notes): add Phase 9 User Guide and Developer Note notes`

## Group 7 — Bookkeeping

- Bump `CMakeLists.txt` `project(gridcalc VERSION 0.8.0)` → `0.9.0`.
- Update `test/version_test.cpp` (`ReturnsZeroEightZero` → `ReturnsZeroNineZero`).
- `CHANGELOG.md`: new `## 0.9.0 — Phase 9: Spectral verification harness` entry listing the vendored dependency, the new namespace, the build flag, and the cross-check coverage.
- `specs/STATUS.md`: bump "Last updated"; move Phase 9 row to Done; update "Where the project stands" + "Next action" (Phase 10 — Documentation infrastructure). Add new "Decisions worth knowing" entries: vendored PocketFFT path; hybrid spectral API; verification-only role of `spectral::*`; Nyquist-zeroing convention.

**Commit message:** `chore: bump version to 0.9.0 and refresh STATUS for Phase 9`
