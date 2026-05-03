/// \file
/// \brief Thin wrapper over PocketFFT for forward / backward 3D
///        real-to-complex FFTs on `gridcalc::core::Field<double>`.
///
/// **Verification only.** The spectral path is not the production
/// differentiation engine — it exists in `gridcalc::spectral` as an
/// independent reference against which the periodic FD operators are
/// cross-checked in CI. See `specs/mission.md` and `specs/roadmap.md`
/// Phase 9.
///
/// Layout: matches Phase 1's locked-in i-fastest layout for `Field<double>`
/// (`linear_index = i + Nx * (j + Ny * k)`). The complex spectrum returned
/// by `forwardR2C` has shape `(Nx/2 + 1) × Ny × Nz` with the same i-fastest
/// indexing in the reduced first dimension. PocketFFT is invoked with
/// column-major (i-fastest) byte strides and `axes = {2, 1, 0}` so that
/// the rfft compresses axis 0 (the x-axis).
/// \since 0.9.0

#pragma once

#ifdef GRIDCALC_USE_FFT

#include <complex>
#include <cstddef>
#include <vector>

#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>

#include <pocketfft_hdronly.h>

namespace gridcalc::spectral {

/// \brief Complex 3D spectrum container, i-fastest layout, half-compressed
///        on the x-axis (rfft).
///
/// Size is `(Nx/2 + 1) * Ny * Nz`. Linear index of element `(nx, ny, nz)`:
/// `nx + (Nx/2 + 1) * (ny + Ny * nz)`.
/// \since 0.9.0
using ComplexField = std::vector<std::complex<double>>;

/// \brief Returns the size of the half-compressed x-axis: `Nx / 2 + 1`.
/// \param grid  Source grid.
/// \returns The compressed first-axis dimension.
/// \since 0.9.0
inline std::size_t rfftXSize(const core::Grid& grid) {
  return static_cast<std::size_t>(grid.getNx() / 2 + 1);
}

/// \brief Returns the linear index `nx + (Nx/2+1) * (ny + Ny * nz)` for a
///        spectrum element.
/// \param grid  Source grid (provides `Nx`, `Ny`).
/// \param nx    X-axis spectrum index in `[0, Nx/2]`.
/// \param ny    Y-axis spectrum index in `[0, Ny)`.
/// \param nz    Z-axis spectrum index in `[0, Nz)`.
/// \returns Linear index into a `ComplexField`.
/// \since 0.9.0
inline std::size_t spectrumIndex(const core::Grid& grid, std::size_t nx,
                                 std::size_t ny, std::size_t nz) {
  const std::size_t nx_size = rfftXSize(grid);
  const std::size_t Ny = static_cast<std::size_t>(grid.getNy());
  return nx + nx_size * (ny + Ny * nz);
}

/// \brief Forward 3D real-to-complex FFT on a `Field<double>`.
///
/// Unnormalized: the standard convention is to divide by $N_x N_y N_z$
/// only on the inverse transform (`backwardC2R`). PocketFFT is invoked
/// with column-major byte strides matching `Field`'s i-fastest layout.
/// \param field  Real input field.
/// \returns A `ComplexField` of size `(Nx/2 + 1) * Ny * Nz`, i-fastest layout.
/// \since 0.9.0
inline ComplexField forwardR2C(const core::Field<double>& field) {
  const auto& grid = field.getGrid();
  const std::size_t Nx = static_cast<std::size_t>(grid.getNx());
  const std::size_t Ny = static_cast<std::size_t>(grid.getNy());
  const std::size_t Nz = static_cast<std::size_t>(grid.getNz());
  const std::size_t nx_out = rfftXSize(grid);

  const pocketfft::shape_t shape_in = {Nx, Ny, Nz};
  const pocketfft::stride_t stride_in = {
      static_cast<std::ptrdiff_t>(sizeof(double)),
      static_cast<std::ptrdiff_t>(Nx * sizeof(double)),
      static_cast<std::ptrdiff_t>(Nx * Ny * sizeof(double))};
  const pocketfft::stride_t stride_out = {
      static_cast<std::ptrdiff_t>(sizeof(std::complex<double>)),
      static_cast<std::ptrdiff_t>(nx_out * sizeof(std::complex<double>)),
      static_cast<std::ptrdiff_t>(nx_out * Ny * sizeof(std::complex<double>))};
  const pocketfft::shape_t axes = {2, 1, 0};

  ComplexField spectrum(nx_out * Ny * Nz);
  pocketfft::r2c(shape_in, stride_in, stride_out, axes, /*forward=*/true,
                 field.data(), spectrum.data(), /*fct=*/1.0);
  return spectrum;
}

/// \brief Backward 3D complex-to-real FFT on a `ComplexField`.
///
/// Normalized: divides the result by $N_x N_y N_z$ so that
/// `backwardC2R(forwardR2C(field))` equals `field` to round-off.
/// \param spectrum  Half-compressed complex spectrum (i-fastest, shape
///                  `(Nx/2 + 1) × Ny × Nz`).
/// \param grid      Target grid; provides $N_x, N_y, N_z, h_a$.
/// \returns A real `Field<double>` on the supplied grid.
/// \since 0.9.0
inline core::Field<double> backwardC2R(const ComplexField& spectrum,
                                       const core::Grid& grid) {
  const std::size_t Nx = static_cast<std::size_t>(grid.getNx());
  const std::size_t Ny = static_cast<std::size_t>(grid.getNy());
  const std::size_t Nz = static_cast<std::size_t>(grid.getNz());
  const std::size_t nx_in = rfftXSize(grid);
  const double inv_total = 1.0 / static_cast<double>(Nx * Ny * Nz);

  const pocketfft::shape_t shape_out = {Nx, Ny, Nz};
  const pocketfft::stride_t stride_in = {
      static_cast<std::ptrdiff_t>(sizeof(std::complex<double>)),
      static_cast<std::ptrdiff_t>(nx_in * sizeof(std::complex<double>)),
      static_cast<std::ptrdiff_t>(nx_in * Ny * sizeof(std::complex<double>))};
  const pocketfft::stride_t stride_out = {
      static_cast<std::ptrdiff_t>(sizeof(double)),
      static_cast<std::ptrdiff_t>(Nx * sizeof(double)),
      static_cast<std::ptrdiff_t>(Nx * Ny * sizeof(double))};
  const pocketfft::shape_t axes = {2, 1, 0};

  core::Field<double> result(grid);
  pocketfft::c2r(shape_out, stride_in, stride_out, axes, /*forward=*/false,
                 spectrum.data(), result.data(), /*fct=*/inv_total);
  return result;
}

}  // namespace gridcalc::spectral

#endif  // GRIDCALC_USE_FFT
