/// \file
/// \brief Declares `gridcalc::core::Grid`, a Cartesian-orthogonal periodic
///        sampling grid in 3D.
/// \since 0.1.0

#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>

#include <gridcalc/core/EigenAliases.hpp>

namespace gridcalc::core {

/// \brief A Cartesian-orthogonal sampling grid: integer extents
///        `(Nx, Ny, Nz)` and a per-axis cell size `(hx, hy, hz)`.
///
/// Sample point `(i, j, k)` sits at Cartesian position
/// `(i * hx, j * hy, k * hz)`; the grid origin is at `(0, 0, 0)` (no half-cell
/// offset). The total domain extent is `(Nx*hx, Ny*hy, Nz*hz)`. Storage layout
/// for any `Field` built on this Grid is i-fastest:
/// `linear_index = i + Nx * (j + Ny * k)`.
///
/// Non-orthogonal Bravais lattices arrive in Phase 15; multi-atom basis in
/// Phase 16. Phase 1 keeps this strictly Cartesian.
/// \since 0.1.0
class Grid {
 public:
  /// \brief Constructs a Grid from extents and per-axis cell sizes.
  /// \param Nx         Number of sample points along x. Must be `> 0`.
  /// \param Ny         Number of sample points along y. Must be `> 0`.
  /// \param Nz         Number of sample points along z. Must be `> 0`.
  /// \param cell_size  Per-axis cell sizes `(hx, hy, hz)`. Every component
  ///                   must be strictly positive.
  /// \throws std::invalid_argument if any extent is non-positive or any
  ///         cell-size component is non-positive.
  /// \since 0.1.0
  Grid(int Nx, int Ny, int Nz, const Vec3d& cell_size)
      : _Nx(Nx), _Ny(Ny), _Nz(Nz), _cell_size(cell_size) {
    if (Nx <= 0 || Ny <= 0 || Nz <= 0) {
      throw std::invalid_argument("Grid: extents (Nx, Ny, Nz) must be > 0; got (" +
                                  std::to_string(Nx) + ", " + std::to_string(Ny) + ", " +
                                  std::to_string(Nz) + ")");
    }
    if (cell_size(0) <= 0.0 || cell_size(1) <= 0.0 || cell_size(2) <= 0.0) {
      throw std::invalid_argument("Grid: cell_size components must be > 0; got (" +
                                  std::to_string(cell_size(0)) + ", " +
                                  std::to_string(cell_size(1)) + ", " +
                                  std::to_string(cell_size(2)) + ")");
    }
  }

  /// \brief Returns the number of sample points along x.
  /// \returns `Nx`.
  /// \since 0.1.0
  int getNx() const noexcept { return _Nx; }

  /// \brief Returns the number of sample points along y.
  /// \returns `Ny`.
  /// \since 0.1.0
  int getNy() const noexcept { return _Ny; }

  /// \brief Returns the number of sample points along z.
  /// \returns `Nz`.
  /// \since 0.1.0
  int getNz() const noexcept { return _Nz; }

  /// \brief Returns the per-axis cell sizes `(hx, hy, hz)`.
  /// \returns Const reference to the stored `Vec3d`.
  /// \since 0.1.0
  const Vec3d& getCellSize() const noexcept { return _cell_size; }

  /// \brief Returns the cell volume `hx * hy * hz`.
  /// \returns Cartesian volume of one grid cell.
  /// \since 0.1.0
  double getCellVolume() const noexcept {
    return _cell_size(0) * _cell_size(1) * _cell_size(2);
  }

  /// \brief Returns the total number of sample points `Nx * Ny * Nz`.
  /// \returns Sample-point count as `std::size_t`.
  /// \since 0.1.0
  std::size_t getSize() const noexcept {
    return static_cast<std::size_t>(_Nx) * static_cast<std::size_t>(_Ny) *
           static_cast<std::size_t>(_Nz);
  }

  /// \brief STL-protocol shim for `getSize()`.
  /// \returns Same as `getSize()`.
  /// \since 0.1.0
  std::size_t size() const noexcept { return getSize(); }

  /// \brief Computes the linear storage offset for a valid `(i, j, k)`.
  ///
  /// Layout is i-fastest: `index = i + Nx * (j + Ny * k)`. Caller must
  /// ensure `i in [0, Nx)`, `j in [0, Ny)`, `k in [0, Nz)`. Wrapping for
  /// out-of-range coordinates is the responsibility of `IndexPolicy`,
  /// applied by `Field::operator()` before this method is invoked.
  /// \param i  In-range x-coordinate.
  /// \param j  In-range y-coordinate.
  /// \param k  In-range z-coordinate.
  /// \returns Linear offset into the i-fastest storage layout.
  /// \since 0.1.0
  std::size_t getLinearIndex(int i, int j, int k) const noexcept {
    return static_cast<std::size_t>(i) +
           static_cast<std::size_t>(_Nx) *
               (static_cast<std::size_t>(j) + static_cast<std::size_t>(_Ny) *
                                                  static_cast<std::size_t>(k));
  }

 private:
  /// \brief Number of sample points along x.
  int _Nx;
  /// \brief Number of sample points along y.
  int _Ny;
  /// \brief Number of sample points along z.
  int _Nz;
  /// \brief Per-axis cell sizes `(hx, hy, hz)`.
  Vec3d _cell_size;
};

}  // namespace gridcalc::core
