/// \file
/// \brief Index-wrapping policies for `gridcalc::core::Field`.
///
/// At Phase 1 only `Periodic` is provided. `Dirichlet` and `Neumann` policies
/// are added in Phase 18.
/// \since 0.1.0

#pragma once

namespace gridcalc::core::IndexPolicy {

/// \brief Periodic-boundary index policy: wraps any integer coordinate into
///        the half-open range `[0, N)` via Euclidean modulo.
///
/// `Periodic` is stateless; the wrap is a pure function of `(i, N)`. The
/// policy is the default template argument of `Field<T, Policy>`.
/// \since 0.1.0
struct Periodic {
  /// \brief Wraps a possibly-out-of-range coordinate into `[0, N)`.
  ///
  /// Implements Euclidean modulo: `wrap(-1, 5) == 4`, `wrap(5, 5) == 0`,
  /// `wrap(-7, 5) == 3`. Required for `N > 0`; behavior is undefined when
  /// `N <= 0` (callers ensure non-empty grids by construction).
  /// \param i  Integer coordinate, possibly negative or `>= N`.
  /// \param N  Axis extent, strictly positive.
  /// \returns The wrapped coordinate in `[0, N)`.
  /// \since 0.1.0
  static constexpr int wrap(int i, int N) noexcept { return ((i % N) + N) % N; }
};

}  // namespace gridcalc::core::IndexPolicy
