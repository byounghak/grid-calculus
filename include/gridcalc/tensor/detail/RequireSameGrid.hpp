/// \file
/// \brief Grid-equality precondition shared by the eager `tensor::*`
///        contraction primitives and the contraction-ET binary nodes.
///
/// `tensor::singleContract`, `tensor::doubleContract`, and the rank-2
/// expression-template binary nodes (`Plus`, `Minus`, `Mul`,
/// `SingleContract`, `DoubleContract`) loop over flat indices `idx`
/// (eager) or evaluate operands at the same `(i, j, k)` (ET) without
/// otherwise checking that both operands sit on the same `Grid`. If
/// they don't, two failure modes appear: equal-`getSize()` /
/// different-extent grids pair physically unrelated cells; smaller
/// right-hand-side storage triggers an out-of-bounds read on the eager
/// path, and the periodic-wrap policy on `Field::operator()` masks the
/// OOB on the ET path while still pairing wrong cells. This helper
/// surfaces both failure modes loudly with a descriptive error.
/// \since 0.15.0

#pragma once

#include <stdexcept>
#include <string>

#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>

namespace gridcalc::tensor::detail {

/// \brief Throws `std::invalid_argument` if `lhs.getGrid()` does not
///        bit-match `rhs.getGrid()`.
///
/// Bit-exact comparison via `core::Grid::operator==` (introduced in
/// `0.14.2`). The message names the caller-supplied `context` string
/// and prints both `(Nx, Ny, Nz, hx, hy, hz)` tuples.
///
/// \tparam T        Element type of the operand fields (typically `core::Mat3d`
///                  for the rank-2 contraction primitives, but the helper is
///                  agnostic since the comparison is only on the underlying `Grid`).
/// \param lhs       Left operand; defines the contract.
/// \param rhs       Right operand; must match `lhs.getGrid()`.
/// \param context   Short tag describing the call site (e.g.,
///                  `"tensor::singleContract"`,
///                  `"tensor::expr::DoubleContract"`).
/// \throws std::invalid_argument if `lhs.getGrid() != rhs.getGrid()`.
/// \since 0.15.0
template <class T>
inline void requireSameGrid(const core::Field<T>& lhs,
                            const core::Field<T>& rhs,
                            const char* context) {
    const auto& a = lhs.getGrid();
    const auto& b = rhs.getGrid();
    if (a == b) return;
    const auto& ha = a.getCellSize();
    const auto& hb = b.getCellSize();
    throw std::invalid_argument(
        std::string("gridcalc: ") + context +
        " requires both operands to share the same Grid, but they do not. "
        "lhs (Nx, Ny, Nz, hx, hy, hz) = (" +
        std::to_string(a.getNx()) + ", " + std::to_string(a.getNy()) + ", " +
        std::to_string(a.getNz()) + ", " + std::to_string(ha(0)) + ", " +
        std::to_string(ha(1)) + ", " + std::to_string(ha(2)) +
        "); rhs = (" + std::to_string(b.getNx()) + ", " +
        std::to_string(b.getNy()) + ", " + std::to_string(b.getNz()) + ", " +
        std::to_string(hb(0)) + ", " + std::to_string(hb(1)) + ", " +
        std::to_string(hb(2)) +
        "). The flat-index loop would silently misalign or read past the "
        "end of the smaller buffer.");
}

/// \brief `Grid`-only overload --- throws if two grids do not bit-match.
///
/// Used by the ET combinator nodes (`Plus`, `Minus`, `Mul`, etc.) where
/// the operands' grids may be queried via the node's `grid()` accessor
/// rather than from `Field<T>::getGrid()` directly.
/// \param lhs      Expected grid (typically the lhs operand's grid).
/// \param rhs      Actual grid (typically the rhs operand's grid).
/// \param context  Short tag describing the call site.
/// \throws std::invalid_argument if `lhs != rhs`.
/// \since 0.15.0
inline void requireSameGrid(const core::Grid& lhs,
                            const core::Grid& rhs,
                            const char* context) {
    if (lhs == rhs) return;
    const auto& ha = lhs.getCellSize();
    const auto& hb = rhs.getCellSize();
    throw std::invalid_argument(
        std::string("gridcalc: ") + context +
        " requires both operands to share the same Grid, but they do not. "
        "lhs (Nx, Ny, Nz, hx, hy, hz) = (" +
        std::to_string(lhs.getNx()) + ", " + std::to_string(lhs.getNy()) + ", " +
        std::to_string(lhs.getNz()) + ", " + std::to_string(ha(0)) + ", " +
        std::to_string(ha(1)) + ", " + std::to_string(ha(2)) +
        "); rhs = (" + std::to_string(rhs.getNx()) + ", " +
        std::to_string(rhs.getNy()) + ", " + std::to_string(rhs.getNz()) + ", " +
        std::to_string(hb(0)) + ", " + std::to_string(hb(1)) + ", " +
        std::to_string(hb(2)) +
        ").");
}

}  // namespace gridcalc::tensor::detail
