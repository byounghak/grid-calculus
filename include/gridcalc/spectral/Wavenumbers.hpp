/// \file
/// \brief Per-axis wavenumber arrays for FFT-based spectral derivatives.
///
/// **Verification only.** The spectral path is not the production
/// differentiation engine — it exists in `gridcalc::spectral` as an
/// independent reference against which the periodic FD operators are
/// cross-checked in CI. See `specs/mission.md` and `specs/roadmap.md`
/// Phase 9.
///
/// Wavenumber sign convention is the standard `numpy.fft`-style: for a
/// full spectrum of length `N` the wavenumber at index `n` is
/// $k[n] = 2\pi\,n^{\text{signed}} / (N\,h)$ where
/// $n^{\text{signed}} = n$ if $2n < N$ else $n - N$.
///
/// The boundary `2n < N` (rather than `n < N/2`) is required so the
/// convention applies uniformly to even and odd `N`:
///  - **Even `N`**: indices `[0, ..., N/2 - 1]` carry the positive
///    harmonics; index `N/2` is the Nyquist mode and lands in the
///    negative branch with $n^{\text{signed}} = -N/2$.
///  - **Odd `N`**: there is no Nyquist mode. Indices
///    `[0, ..., (N - 1)/2]` carry the positive harmonics (highest
///    is `(N - 1)/2`); indices `[(N + 1)/2, ..., N - 1]` carry the
///    negative harmonics.
///
/// The rfft (x-axis) half-spectrum carries only non-negative
/// frequencies $n \in [0, N/2]$ (length `N/2 + 1`); the sign-flip
/// path is absent and the rfft form is correct for both parities.
/// \since 0.9.0

#pragma once

#ifdef GRIDCALC_USE_FFT

#include <cstddef>
#include <vector>

#include <gridcalc/core/Grid.hpp>

namespace gridcalc::spectral {

/// \brief Returns the rfft half-spectrum wavenumber array along the x-axis.
///
/// Shape: `Nx/2 + 1`. Entries are
/// $k_x[n_x] = 2\pi n_x / (N_x h_x)$ for $n_x \in [0, N_x/2]$.
/// Used by the spectral derivative operators that consume the half-spectrum
/// produced by `forwardR2C`.
/// \param grid  Source grid; provides $N_x$ and $h_x$.
/// \returns A `std::vector<double>` of length `Nx/2 + 1`.
/// \since 0.9.0
inline std::vector<double> kxRfft(const core::Grid& grid) {
  constexpr double kTwoPi = 6.283185307179586476925286766559;
  const int Nx = grid.getNx();
  const double hx = grid.getCellSize()(0);
  const double scale = kTwoPi / (static_cast<double>(Nx) * hx);
  std::vector<double> kx(static_cast<std::size_t>(Nx / 2 + 1));
  for (std::size_t n = 0; n < kx.size(); ++n) {
    kx[n] = scale * static_cast<double>(n);
  }
  return kx;
}

/// \brief Returns the full-spectrum wavenumber array along the y-axis.
///
/// Shape: `Ny`. Entries follow the standard signed `numpy.fft.fftfreq`
/// convention: $k_y[n_y] = 2\pi n_y^{\text{signed}} / (N_y h_y)$
/// where $n_y^{\text{signed}} = n_y$ if $2 n_y < N_y$ else $n_y - N_y$.
/// The boundary form `2n < N` (rather than `n < N/2`) handles even
/// and odd `Ny` uniformly — see the file-level block.
/// \param grid  Source grid; provides $N_y$ and $h_y$.
/// \returns A `std::vector<double>` of length `Ny`.
/// \since 0.9.0 (function); 0.14.3 (odd-`N` off-by-one fix).
inline std::vector<double> kyFull(const core::Grid& grid) {
  constexpr double kTwoPi = 6.283185307179586476925286766559;
  const int Ny = grid.getNy();
  const double hy = grid.getCellSize()(1);
  const double scale = kTwoPi / (static_cast<double>(Ny) * hy);
  std::vector<double> ky(static_cast<std::size_t>(Ny));
  for (int n = 0; n < Ny; ++n) {
    const int n_signed = (2 * n < Ny) ? n : n - Ny;
    ky[static_cast<std::size_t>(n)] = scale * static_cast<double>(n_signed);
  }
  return ky;
}

/// \brief Returns the full-spectrum wavenumber array along the z-axis.
///
/// Shape: `Nz`. Same signed `numpy.fft.fftfreq` convention as
/// `kyFull`: $k_z[n_z] = 2\pi n_z^{\text{signed}} / (N_z h_z)$ where
/// $n_z^{\text{signed}} = n_z$ if $2 n_z < N_z$ else $n_z - N_z$. The
/// boundary form `2n < N` handles even and odd `Nz` uniformly — see
/// the file-level block.
/// \param grid  Source grid; provides $N_z$ and $h_z$.
/// \returns A `std::vector<double>` of length `Nz`.
/// \since 0.9.0 (function); 0.14.3 (odd-`N` off-by-one fix).
inline std::vector<double> kzFull(const core::Grid& grid) {
  constexpr double kTwoPi = 6.283185307179586476925286766559;
  const int Nz = grid.getNz();
  const double hz = grid.getCellSize()(2);
  const double scale = kTwoPi / (static_cast<double>(Nz) * hz);
  std::vector<double> kz(static_cast<std::size_t>(Nz));
  for (int n = 0; n < Nz; ++n) {
    const int n_signed = (2 * n < Nz) ? n : n - Nz;
    kz[static_cast<std::size_t>(n)] = scale * static_cast<double>(n_signed);
  }
  return kz;
}

}  // namespace gridcalc::spectral

#endif  // GRIDCALC_USE_FFT
