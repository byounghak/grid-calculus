/// \file
/// \brief FFT-based spectral derivative operators on periodic grids.
///
/// **Verification only.** The spectral path is not the production
/// differentiation engine — it exists in `gridcalc::spectral` as an
/// independent reference against which the periodic FD operators in
/// `gridcalc::diff` are cross-checked in CI. See `specs/mission.md`
/// and `specs/roadmap.md` Phase 9.
///
/// Provides:
/// - `spectral::partial<Nx, Ny, Nz>(field)` — generic multi-index spectral
///   partial. Multiplies the spectrum by $(ik_x)^{N_x}(ik_y)^{N_y}(ik_z)^{N_z}$
///   and inverse-transforms.
/// - `spectral::laplacian(field)` — fused $-|\mathbf{k}|^2$ multiplication.
/// - `spectral::biharmonic(field)` — fused $|\mathbf{k}|^4$ multiplication.
/// - `spectral::gradient(field)` — three first-derivatives packed into
///   `Field<Vec3d>`; one forward + three backward FFTs.
///
/// Nyquist convention: for any axis with an **odd** per-axis derivative
/// count **and an even axis extent `N_a`**, the mode at
/// $n_a = N_a / 2$ (the Nyquist mode for even `N_a`) is set to zero
/// before applying the spectral multiplier. This is the standard remedy
/// for the $\pm\pi/h$ sign ambiguity at the Nyquist that would
/// otherwise alias under the inverse transform.
///
/// **Odd `N_a` has no Nyquist mode by construction** — the highest
/// positive harmonic is at $(N_a - 1)/2$ and is a regular real
/// harmonic. The zeroing path is short-circuited by the parity gate
/// (`N_a % 2 == 0`); the highest positive harmonic is preserved with
/// its full multiplier. Verified by `test/spectral_partial_nyquist_test.cpp`.
/// \since 0.9.0

#pragma once

#ifdef GRIDCALC_USE_FFT

#include <complex>
#include <cstddef>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/spectral/Fft.hpp>
#include <gridcalc/spectral/Wavenumbers.hpp>

namespace gridcalc::spectral {

namespace detail {

/// \brief `constexpr` non-negative integer power: `v^n` for `n >= 0`.
/// `0^0` is defined to be `1`.
/// \since 0.9.0
constexpr double powInt(double v, int n) {
  double r = 1.0;
  for (int i = 0; i < n; ++i) r *= v;
  return r;
}

}  // namespace detail

/// \brief Returns the multi-index partial derivative
///        $\partial^{N_x+N_y+N_z}\psi / \partial x^{N_x}\partial y^{N_y}\partial z^{N_z}$
///        via spectral multiplication by $(ik_x)^{N_x}(ik_y)^{N_y}(ik_z)^{N_z}$.
///
/// For each spectrum element the multiplier factors as
/// $i^{N_x+N_y+N_z}\,k_x^{N_x}\,k_y^{N_y}\,k_z^{N_z}$. The leading
/// $i^{N}$ is real or purely imaginary depending on $N \bmod 4$.
/// On any axis with an **odd** derivative count **and an even axis
/// extent `N_a`**, the mode at $n_a = N_a / 2$ (the Nyquist) is
/// zeroed. Odd `N_a` has no Nyquist mode and the highest positive
/// harmonic at $(N_a - 1)/2$ keeps its full multiplier.
/// \tparam Nx     X-axis derivative count ($\geq 0$).
/// \tparam Ny     Y-axis derivative count ($\geq 0$).
/// \tparam Nz     Z-axis derivative count ($\geq 0$).
/// \param field   Real input field.
/// \returns A new `Field<double>` holding the requested partial.
/// \since 0.9.0
template <int Nx, int Ny, int Nz>
inline core::Field<double> partial(const core::Field<double>& field) {
  const auto& grid = field.getGrid();
  ComplexField spectrum = forwardR2C(field);

  const auto kx = kxRfft(grid);
  const auto ky = kyFull(grid);
  const auto kz = kzFull(grid);

  const std::size_t nx_size = rfftXSize(grid);
  const std::size_t Ny_size = static_cast<std::size_t>(grid.getNy());
  const std::size_t Nz_size = static_cast<std::size_t>(grid.getNz());
  const int Nx_grid = grid.getNx();
  const int Ny_grid = grid.getNy();
  const int Nz_grid = grid.getNz();
  constexpr int total = Nx + Ny + Nz;

  std::complex<double> i_power;
  switch (total % 4) {
    case 0: i_power = {1.0, 0.0}; break;
    case 1: i_power = {0.0, 1.0}; break;
    case 2: i_power = {-1.0, 0.0}; break;
    case 3: i_power = {0.0, -1.0}; break;
    default: i_power = {1.0, 0.0};
  }

  for (std::size_t nz = 0; nz < Nz_size; ++nz) {
    const double kz_pow = detail::powInt(kz[nz], Nz);
    for (std::size_t ny = 0; ny < Ny_size; ++ny) {
      const double ky_pow = detail::powInt(ky[ny], Ny);
      for (std::size_t nx = 0; nx < nx_size; ++nx) {
        const double kx_pow = detail::powInt(kx[nx], Nx);
        std::complex<double> mult = i_power * (kx_pow * ky_pow * kz_pow);

        if constexpr (Nx % 2 == 1) {
          if (Nx_grid % 2 == 0 && nx == nx_size - 1) mult = {0.0, 0.0};
        }
        if constexpr (Ny % 2 == 1) {
          if (Ny_grid % 2 == 0 &&
              ny == static_cast<std::size_t>(Ny_grid / 2)) {
            mult = {0.0, 0.0};
          }
        }
        if constexpr (Nz % 2 == 1) {
          if (Nz_grid % 2 == 0 &&
              nz == static_cast<std::size_t>(Nz_grid / 2)) {
            mult = {0.0, 0.0};
          }
        }

        spectrum[spectrumIndex(grid, nx, ny, nz)] *= mult;
      }
    }
  }

  return backwardC2R(spectrum, grid);
}

/// \brief Returns $\nabla^2\psi$ via fused spectral multiplication by
///        $-|\mathbf{k}|^2$.
///
/// Single forward + single backward FFT (vs three of each for the
/// composition `partial<2,0,0> + partial<0,2,0> + partial<0,0,2>`).
/// \param field  Real input field.
/// \returns A new `Field<double>` holding $\nabla^2\psi$.
/// \since 0.9.0
inline core::Field<double> laplacian(const core::Field<double>& field) {
  const auto& grid = field.getGrid();
  ComplexField spectrum = forwardR2C(field);
  const auto kx = kxRfft(grid);
  const auto ky = kyFull(grid);
  const auto kz = kzFull(grid);
  const std::size_t nx_size = rfftXSize(grid);
  const std::size_t Ny_size = static_cast<std::size_t>(grid.getNy());
  const std::size_t Nz_size = static_cast<std::size_t>(grid.getNz());

  for (std::size_t nz = 0; nz < Nz_size; ++nz) {
    const double kz_sq = kz[nz] * kz[nz];
    for (std::size_t ny = 0; ny < Ny_size; ++ny) {
      const double ky_sq = ky[ny] * ky[ny];
      for (std::size_t nx = 0; nx < nx_size; ++nx) {
        const double kx_sq = kx[nx] * kx[nx];
        const double mult = -(kx_sq + ky_sq + kz_sq);
        spectrum[spectrumIndex(grid, nx, ny, nz)] *= mult;
      }
    }
  }
  return backwardC2R(spectrum, grid);
}

/// \brief Returns $\nabla^4\psi$ via fused spectral multiplication by
///        $|\mathbf{k}|^4$.
///
/// Single forward + single backward FFT.
/// \param field  Real input field.
/// \returns A new `Field<double>` holding $\nabla^4\psi$.
/// \since 0.9.0
inline core::Field<double> biharmonic(const core::Field<double>& field) {
  const auto& grid = field.getGrid();
  ComplexField spectrum = forwardR2C(field);
  const auto kx = kxRfft(grid);
  const auto ky = kyFull(grid);
  const auto kz = kzFull(grid);
  const std::size_t nx_size = rfftXSize(grid);
  const std::size_t Ny_size = static_cast<std::size_t>(grid.getNy());
  const std::size_t Nz_size = static_cast<std::size_t>(grid.getNz());

  for (std::size_t nz = 0; nz < Nz_size; ++nz) {
    const double kz_sq = kz[nz] * kz[nz];
    for (std::size_t ny = 0; ny < Ny_size; ++ny) {
      const double ky_sq = ky[ny] * ky[ny];
      for (std::size_t nx = 0; nx < nx_size; ++nx) {
        const double kx_sq = kx[nx] * kx[nx];
        const double k_sq = kx_sq + ky_sq + kz_sq;
        spectrum[spectrumIndex(grid, nx, ny, nz)] *= (k_sq * k_sq);
      }
    }
  }
  return backwardC2R(spectrum, grid);
}

/// \brief Returns the spectral gradient $\nabla\psi$ as a `Field<Vec3d>`.
///
/// One forward FFT plus three backward FFTs (one per Cartesian component).
/// Each component applies its own $i k_a$ spectrum multiplier with
/// per-axis Nyquist zeroing.
/// \param field  Real scalar input field.
/// \returns A new `Field<Vec3d>` holding $\nabla\psi$.
/// \since 0.9.0
inline core::Field<core::Vec3d> gradient(const core::Field<double>& field) {
  const auto& grid = field.getGrid();
  ComplexField base = forwardR2C(field);
  const auto kx = kxRfft(grid);
  const auto ky = kyFull(grid);
  const auto kz = kzFull(grid);
  const std::size_t nx_size = rfftXSize(grid);
  const std::size_t Ny_size = static_cast<std::size_t>(grid.getNy());
  const std::size_t Nz_size = static_cast<std::size_t>(grid.getNz());
  const int Nx_grid = grid.getNx();
  const int Ny_grid = grid.getNy();
  const int Nz_grid = grid.getNz();

  ComplexField spec_x = base;
  ComplexField spec_y = base;
  ComplexField spec_z = base;
  for (std::size_t nz = 0; nz < Nz_size; ++nz) {
    for (std::size_t ny = 0; ny < Ny_size; ++ny) {
      for (std::size_t nx = 0; nx < nx_size; ++nx) {
        const std::size_t idx = spectrumIndex(grid, nx, ny, nz);
        std::complex<double> mx{0.0, kx[nx]};
        std::complex<double> my{0.0, ky[ny]};
        std::complex<double> mz{0.0, kz[nz]};
        if (Nx_grid % 2 == 0 && nx == nx_size - 1) mx = {0.0, 0.0};
        if (Ny_grid % 2 == 0 &&
            ny == static_cast<std::size_t>(Ny_grid / 2)) {
          my = {0.0, 0.0};
        }
        if (Nz_grid % 2 == 0 &&
            nz == static_cast<std::size_t>(Nz_grid / 2)) {
          mz = {0.0, 0.0};
        }
        spec_x[idx] *= mx;
        spec_y[idx] *= my;
        spec_z[idx] *= mz;
      }
    }
  }

  core::Field<double> dx_field = backwardC2R(spec_x, grid);
  core::Field<double> dy_field = backwardC2R(spec_y, grid);
  core::Field<double> dz_field = backwardC2R(spec_z, grid);

  core::Field<core::Vec3d> result(grid);
  const std::size_t total = result.getSize();
  for (std::size_t n = 0; n < total; ++n) {
    result.data()[n] = core::Vec3d(dx_field.data()[n], dy_field.data()[n],
                                   dz_field.data()[n]);
  }
  return result;
}

}  // namespace gridcalc::spectral

#endif  // GRIDCALC_USE_FFT
