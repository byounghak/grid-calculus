/// \file
/// \brief Conservative cell-flux discretization of the heterogeneous-D
///        diffusion operator: $\nabla\cdot(D\,\nabla\psi)$.
///
/// At each cell `(i, j, k)`, the operator computes the net flux out of
/// the cell across its six faces, divided by the cell volume. Per axis,
/// the face flux is the standard FVM 2nd-order central-difference:
/// \verbatim
///   D_face_E (psi(i+1) - psi(i)) - D_face_W (psi(i) - psi(i-1))
///   ---------------------------------------------------------------
///                            h_x^2
/// \endverbatim
/// with face-centered diffusivity `D_face` taken as the **harmonic mean**
/// of the two adjacent cell values:
/// \verbatim
///   D_face = 2 * D_left * D_right / (D_left + D_right)
/// \endverbatim
/// The harmonic mean is the textbook FVM choice for diffusion (LeVeque
/// 2007 §4.2; Patankar 1980 §4.2): it recovers the correct steady-state
/// flux when `D` jumps discontinuously between cells, where the
/// arithmetic mean would bias toward the higher-D side.
///
/// **Boundaries.** Periodic only at Phase 14. Non-periodic boundary
/// conditions arrive in Phase 19 and consume the FVM foundation.
/// \since 0.14.0

#pragma once

#include <cstddef>

#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>

namespace gridcalc::fvm {

namespace detail {

/// \brief Harmonic mean of two strictly positive values: $2\,a\,b / (a+b)$.
///
/// Caller ensures `a > 0` and `b > 0`. The harmonic mean is undefined
/// when either is zero; the public `cellLaplacian` contract requires
/// `D > 0` everywhere on the grid.
/// \since 0.14.0
inline double harmonicMean(double a, double b) noexcept {
  return 2.0 * a * b / (a + b);
}

}  // namespace detail

/// \brief Returns $\nabla\cdot(D\,\nabla\psi)$ as a fresh, same-shape
///        `Field<double>`, computed via cell-flux discretization with
///        face-centered harmonic-mean D-averaging.
///
/// Each Cartesian axis contributes
/// $(D_{a,+} (\psi_{a,+}-\psi) - D_{a,-} (\psi-\psi_{a,-}))/h_a^2$
/// where $D_{a,\pm}$ are the harmonic means of `D` across the two faces
/// adjacent to the cell along axis $a$. Wrapping at the domain boundary
/// is the responsibility of the input `Field`'s `Policy` (default
/// `IndexPolicy::Periodic`); `cellLaplacian` itself makes no boundary
/// decisions.
///
/// On a uniform Cartesian grid with constant `D = c`, the harmonic mean
/// reduces to `c` and the output equals `c * diff::laplacian(psi)` to
/// round-off — verified by the agreement gate
/// `fvm::cellLaplacian == diff::laplacian` test in
/// `test/fvm_cell_laplacian_test.cpp`.
///
/// \param psi  Input scalar field.
/// \param D    Diffusivity field on the same `Grid` as `psi`. **Must be
///             strictly positive** at every cell; the harmonic mean is
///             undefined otherwise. Behavior is undefined if either
///             contract is violated.
/// \returns A new `Field<double>` holding $\nabla\cdot(D\,\nabla\psi)$.
/// \since 0.14.0
inline core::Field<double> cellLaplacian(const core::Field<double>& psi,
                                         const core::Field<double>& D) {
  const auto& grid = psi.getGrid();
  const auto& h = grid.getCellSize();
  const double inv_hx2 = 1.0 / (h(0) * h(0));
  const double inv_hy2 = 1.0 / (h(1) * h(1));
  const double inv_hz2 = 1.0 / (h(2) * h(2));

  core::Field<double> result(grid);
  for (int k = 0; k < grid.getNz(); ++k) {
    for (int j = 0; j < grid.getNy(); ++j) {
      for (int i = 0; i < grid.getNx(); ++i) {
        const double psi_c = psi(i, j, k);
        const double D_c = D(i, j, k);

        const double D_face_xp = detail::harmonicMean(D_c, D(i + 1, j, k));
        const double D_face_xm = detail::harmonicMean(D(i - 1, j, k), D_c);
        const double D_face_yp = detail::harmonicMean(D_c, D(i, j + 1, k));
        const double D_face_ym = detail::harmonicMean(D(i, j - 1, k), D_c);
        const double D_face_zp = detail::harmonicMean(D_c, D(i, j, k + 1));
        const double D_face_zm = detail::harmonicMean(D(i, j, k - 1), D_c);

        const double rhs_x = (D_face_xp * (psi(i + 1, j, k) - psi_c) -
                              D_face_xm * (psi_c - psi(i - 1, j, k))) *
                             inv_hx2;
        const double rhs_y = (D_face_yp * (psi(i, j + 1, k) - psi_c) -
                              D_face_ym * (psi_c - psi(i, j - 1, k))) *
                             inv_hy2;
        const double rhs_z = (D_face_zp * (psi(i, j, k + 1) - psi_c) -
                              D_face_zm * (psi_c - psi(i, j, k - 1))) *
                             inv_hz2;
        result(i, j, k) = rhs_x + rhs_y + rhs_z;
      }
    }
  }
  return result;
}

}  // namespace gridcalc::fvm
