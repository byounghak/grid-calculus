/// \file
/// \brief Grid-equality precondition shared by `solve::integrate(...)`
///        and `detail::axpyFresh`.
///
/// `solve::integrate` calls a user-supplied `rhs(psi)` callable each
/// stage and then loops `n = 0..psi.getSize()-1` over the returned
/// field's `data()` pointer. The static-assert
/// `std::is_invocable_r_v<Field<double>, Rhs&, const Field<double>&>`
/// only constrains the *return type*, not the runtime `Grid` carried
/// by the returned field. If the RHS allocates against a different
/// `Grid`, the integrator's flat-index loop silently mixes derivatives
/// from misaligned cells (or reads past the end of the smaller field's
/// storage). This helper lets each stage's entry point reject those
/// returns loudly instead.
/// \since 0.14.2

#pragma once

#include <stdexcept>
#include <string>

#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>

namespace gridcalc::solve::detail {

/// \brief Throws `std::invalid_argument` if `actual.getGrid()` does
///        not bit-match `expected.getGrid()`.
///
/// Bit-exact comparison via `Grid::operator==`. The message names the
/// caller-supplied `context` string (so the user can locate the
/// offending call site without reading source) and prints both the
/// expected and actual `(Nx, Ny, Nz, hx, hy, hz)` tuples.
///
/// \param expected  Field whose `Grid` defines the contract.
/// \param actual    Field whose `Grid` must match.
/// \param context   Short tag describing the call site (e.g.,
///                  `"solve::integrate(... RK4): k2"`,
///                  `"detail::axpyFresh"`).
/// \throws std::invalid_argument if `actual.getGrid() != expected.getGrid()`.
/// \since 0.14.2
inline void requireSameGrid(const core::Field<double>& expected,
                            const core::Field<double>& actual,
                            const char* context) {
  const auto& a = expected.getGrid();
  const auto& b = actual.getGrid();
  if (a == b) return;
  const auto& ha = a.getCellSize();
  const auto& hb = b.getCellSize();
  throw std::invalid_argument(
      std::string("gridcalc: ") + context +
      " produced a Field on a different Grid than expected. "
      "Expected (Nx, Ny, Nz, hx, hy, hz) = (" +
      std::to_string(a.getNx()) + ", " + std::to_string(a.getNy()) + ", " +
      std::to_string(a.getNz()) + ", " + std::to_string(ha(0)) + ", " +
      std::to_string(ha(1)) + ", " + std::to_string(ha(2)) +
      "); got (" + std::to_string(b.getNx()) + ", " +
      std::to_string(b.getNy()) + ", " + std::to_string(b.getNz()) + ", " +
      std::to_string(hb(0)) + ", " + std::to_string(hb(1)) + ", " +
      std::to_string(hb(2)) +
      "). The integrator loops over psi.getSize() against the returned "
      "field's data() pointer; a Grid mismatch silently misaligns "
      "(or reads past the end of) the buffer.");
}

}  // namespace gridcalc::solve::detail
