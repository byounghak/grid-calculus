# User Guide — Phase 9: Spectral verification harness

> Spec dir: `specs/2026-05-03-phase-9-spectral-verification/`. Version: `0.9.0`.

## What this gives you

A second, independent path for computing periodic-grid derivatives via FFT.
**Not a production differentiation engine.** The finite-difference operators
in `gridcalc::diff` (Phases 1, 2, 7, 8) are the production code; the
`gridcalc::spectral` namespace exists to be the gold-standard reference
that gates FD changes in CI. Use the spectral path when you want a
verification check, not when you're solving a real problem.

```cpp
diff::laplacian<4>(psi);        // production differentiation
spectral::laplacian(psi);       // verification reference (FFT-based)
```

The library ships a parameterized GoogleTest fixture
([`test/fd_fft_crosscheck_test.cpp`](../../../test/fd_fft_crosscheck_test.cpp))
that already cross-checks every Phase 1–8 FD operator at every advertised
accuracy order against its spectral counterpart on every CI run. Adding a
new FD operator in a future phase means adding one more line to that
fixture; the gate catches accidental accuracy regressions automatically.

## API

Headers under [`include/gridcalc/spectral/`](../../../include/gridcalc/spectral/):

```cpp
namespace gridcalc::spectral {

// Generic multi-index spectral partial. Multiplies the spectrum by
// (i kx)^Nx (i ky)^Ny (i kz)^Nz; zeros the per-axis Nyquist mode for
// any axis with an odd derivative count.
template <int Nx, int Ny, int Nz>
core::Field<double> partial(const core::Field<double>& field);

// Fused single-FFT-pair operators for the contracted / rank-1 cases.
core::Field<double>      laplacian (const core::Field<double>& field);   // -|k|^2
core::Field<double>      biharmonic(const core::Field<double>& field);   //  |k|^4
core::Field<core::Vec3d> gradient  (const core::Field<double>& field);   // i k

// Lower-level FFT primitives (rarely needed directly).
ComplexField forwardR2C  (const core::Field<double>& field);
core::Field<double>
             backwardC2R (const ComplexField& spectrum, const core::Grid& grid);

// Wavenumber helpers.
std::vector<double> kxRfft(const core::Grid& grid);  // length Nx/2 + 1
std::vector<double> kyFull(const core::Grid& grid);  // length Ny (signed)
std::vector<double> kzFull(const core::Grid& grid);  // length Nz (signed)

}  // namespace gridcalc::spectral
```

The Phase 8 d-prefix partials (`dxx`, `d3dx2dy`, `d4dx2dy2`, etc.) do not
have one-for-one named `spectral::*` counterparts. Use `spectral::partial`
with the matching multi-index instead — for example, the spectral analog
of `diff::d3dx2dy<4>(psi)` is `spectral::partial<2, 1, 0>(psi)`.

## Build flag

```
cmake -B build -DGRIDCALC_USE_FFT=ON   # default: builds the spectral harness
cmake -B build -DGRIDCALC_USE_FFT=OFF  # skips spectral headers and CI gate
```

When `OFF`, the `spectral/` headers are no-ops (`#ifdef GRIDCALC_USE_FFT`
guards) and the FD–FFT cross-check fixture is excluded from the test
binary. The FD code path is unaffected.

## Wavenumber and Nyquist convention

- **x-axis (rfft):** half-spectrum of length `Nx/2 + 1`. Wavenumbers
  $k_x[n_x] = 2\pi n_x / (N_x h_x)$ for $n_x \in [0, N_x/2]$.
- **y, z (full c2c):** signed `numpy.fft`-style.
  $k_a[n_a] = 2\pi n_a^{\text{signed}} / (N_a h_a)$ where
  $n_a^{\text{signed}} = n_a$ for $n_a < N_a/2$ else $n_a - N_a$.
- **Nyquist zeroing.** For any axis whose per-axis derivative count is
  **odd**, the mode at $n_a = N_a/2$ (Nyquist) is set to zero before
  applying the spectral multiplier. This handles the standard $\pm \pi/h$
  sign ambiguity at the Nyquist that would otherwise alias under the
  inverse transform.

## Worked example: a CI-style FD–FFT assertion

Treat the spectral path as the reference and assert that an FD operator
agrees with it within the operator's truncation order:

```cpp
#include <cmath>
#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/diff/Biharmonic.hpp>
#include <gridcalc/spectral/Derivatives.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;

constexpr double kPi = 3.14159265358979323846;

const int N = 64;
const double L = 2.0 * kPi;
const double h = L / N;

Grid grid(N, N, N, Vec3d(h, h, h));
Field<double> psi(grid, [](double x, double y, double z) {
  return std::sin(x + y + z);
});

Field<double> fd  = gridcalc::diff::biharmonic<4>(psi);
Field<double> sp  = gridcalc::spectral::biharmonic(psi);

// Compute relative max-norm discrepancy:
double max_diff = 0.0;
double max_ref  = 0.0;
for (std::size_t n = 0; n < psi.getSize(); ++n) {
  max_diff = std::max(max_diff, std::abs(fd.data()[n] - sp.data()[n]));
  max_ref  = std::max(max_ref,  std::abs(sp.data()[n]));
}
const double rel = max_diff / max_ref;
// Order-4 truncation at N=64 on |k|=sqrt(3): rel ~ h^4 * |k|^2 ~ 3e-4.
// A loose 1e-3 bound is plenty for a CI gate.
assert(rel < 1e-3);
```

The cross-check fixture
[`test/fd_fft_crosscheck_test.cpp`](../../../test/fd_fft_crosscheck_test.cpp)
generalizes this pattern to every Phase 1–8 operator and asserts the
log-log slope on a four-point grid sweep.

## What this does *not* do

- **Production differentiation.** Use `gridcalc::diff::*` instead.
- **Non-orthogonal Bravais grids.** Phase 15.
- **Multi-atom-basis grids.** Phase 16.
- **31 named `spectral::*` d-prefix partials.** Use `spectral::partial<Nx, Ny, Nz>(field)` instead.
- **Spectral solver path.** Phase 19's implicit-scheme spec round will reconsider this; the Phase 9 vendor of PocketFFT is verification-only.

For the wavenumber sign convention, the Nyquist-aliasing argument, and
the design choices behind the hybrid API and the vendored PocketFFT, see
[`docs/developer-note/notes/phase-9-spectral-verification.md`](../../developer-note/notes/phase-9-spectral-verification.md).
