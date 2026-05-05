// Tests for the RHS-grid precondition added at 0.14.2.
//
// solve::integrate(...) (RK4 + ExplicitEuler) and detail::axpyFresh
// now reject any field returned by the user's RHS callable that does
// not carry psi.getGrid() (bit-exact). Without the check, the
// integrator's flat-index loop over psi.getSize() against the
// returned field's data() pointer would silently mix derivatives from
// misaligned cells -- or read past the end of the smaller field's
// std::vector storage.
//
// Coverage:
//  - integrate(... RK4) rejects on three flavors of mismatch:
//      (a) different total cell count (OOB-risk case),
//      (b) same total count, different per-axis extents (silent layout drift),
//      (c) same shape, different cell size (silent 1/h drift).
//  - integrate(... ExplicitEuler) rejects on the same three flavors.
//  - detail::axpyFresh rejects directly when its two inputs differ.
//  - Per-stage label sweep on RK4: a stateful RHS that returns a
//    correct grid for k1..k_{i-1} and a wrong grid for k_i confirms
//    the helper names the offending stage (k1/k2/k3/k4) in its message.
//  - Smoke test: a correctly-grid-preserving RHS (diff::laplacian)
//    integrates without throwing through both integrator overloads.

#include <gtest/gtest.h>

#include <string>

#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/diff/Laplacian.hpp>
#include <gridcalc/solve/ExplicitEuler.hpp>
#include <gridcalc/solve/RK4.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
using gridcalc::solve::ExplicitEuler;
using gridcalc::solve::RK4;
using gridcalc::solve::integrate;

namespace {

Field<double> makeStateField() {
  Grid g(8, 8, 8, Vec3d(0.1, 0.1, 0.1));
  return Field<double>(g, 0.0);
}

// RHS that always returns a field on a different (smaller) grid.
struct DifferentTotalCountRhs {
  Grid wrong_grid{4, 4, 4, Vec3d(0.1, 0.1, 0.1)};
  Field<double> operator()(const Field<double>&) const {
    return Field<double>(wrong_grid, 0.0);
  }
};

// RHS that returns a field with the same total cell count (512) but
// different per-axis extents (4 * 16 * 8 = 512 = 8 * 8 * 8).
struct SameSizeDifferentShapeRhs {
  Grid wrong_grid{4, 16, 8, Vec3d(0.1, 0.1, 0.1)};
  Field<double> operator()(const Field<double>&) const {
    return Field<double>(wrong_grid, 0.0);
  }
};

// RHS that returns a field with the same shape but different cell size
// (silent 1/h-drift case).
struct SameShapeDifferentCellSizeRhs {
  Grid wrong_grid{8, 8, 8, Vec3d(0.2, 0.2, 0.2)};
  Field<double> operator()(const Field<double>&) const {
    return Field<double>(wrong_grid, 0.0);
  }
};

}  // namespace

// =============================================================================
// integrate(... RK4) — three mismatch flavors
// =============================================================================

TEST(IntegratorGridMismatchTest, RK4_RejectsDifferentTotalCount) {
  Field<double> psi = makeStateField();
  EXPECT_THROW(integrate(psi, DifferentTotalCountRhs{}, 0.01, 1, RK4{}),
               std::invalid_argument);
}

TEST(IntegratorGridMismatchTest, RK4_RejectsSameSizeDifferentShape) {
  Field<double> psi = makeStateField();
  EXPECT_THROW(integrate(psi, SameSizeDifferentShapeRhs{}, 0.01, 1, RK4{}),
               std::invalid_argument);
}

TEST(IntegratorGridMismatchTest, RK4_RejectsSameShapeDifferentCellSize) {
  Field<double> psi = makeStateField();
  EXPECT_THROW(integrate(psi, SameShapeDifferentCellSizeRhs{}, 0.01, 1, RK4{}),
               std::invalid_argument);
}

// =============================================================================
// integrate(... ExplicitEuler) — three mismatch flavors
// =============================================================================

TEST(IntegratorGridMismatchTest, ExplicitEuler_RejectsDifferentTotalCount) {
  Field<double> psi = makeStateField();
  EXPECT_THROW(
      integrate(psi, DifferentTotalCountRhs{}, 0.01, 1, ExplicitEuler{}),
      std::invalid_argument);
}

TEST(IntegratorGridMismatchTest, ExplicitEuler_RejectsSameSizeDifferentShape) {
  Field<double> psi = makeStateField();
  EXPECT_THROW(
      integrate(psi, SameSizeDifferentShapeRhs{}, 0.01, 1, ExplicitEuler{}),
      std::invalid_argument);
}

TEST(IntegratorGridMismatchTest, ExplicitEuler_RejectsSameShapeDifferentCellSize) {
  Field<double> psi = makeStateField();
  EXPECT_THROW(
      integrate(psi, SameShapeDifferentCellSizeRhs{}, 0.01, 1, ExplicitEuler{}),
      std::invalid_argument);
}

// =============================================================================
// axpyFresh — direct precondition coverage
// =============================================================================

TEST(IntegratorGridMismatchTest, AxpyFreshRejectsMismatchedShape) {
  Field<double> a(Grid(8, 8, 8, Vec3d(0.1, 0.1, 0.1)), 1.0);
  Field<double> b(Grid(4, 4, 4, Vec3d(0.1, 0.1, 0.1)), 1.0);
  EXPECT_THROW(gridcalc::solve::detail::axpyFresh(a, 0.5, b),
               std::invalid_argument);
}

TEST(IntegratorGridMismatchTest, AxpyFreshRejectsMismatchedCellSize) {
  Field<double> a(Grid(8, 8, 8, Vec3d(0.1, 0.1, 0.1)), 1.0);
  Field<double> b(Grid(8, 8, 8, Vec3d(0.2, 0.1, 0.1)), 1.0);
  EXPECT_THROW(gridcalc::solve::detail::axpyFresh(a, 0.5, b),
               std::invalid_argument);
}

TEST(IntegratorGridMismatchTest, AxpyFreshAcceptsMatchedGrid) {
  Field<double> a(Grid(8, 8, 8, Vec3d(0.1, 0.1, 0.1)), 1.0);
  Field<double> b(Grid(8, 8, 8, Vec3d(0.1, 0.1, 0.1)), 2.0);
  EXPECT_NO_THROW({
    auto out = gridcalc::solve::detail::axpyFresh(a, 0.5, b);
    EXPECT_EQ(out.getSize(), a.getSize());
  });
}

// =============================================================================
// Per-stage label sweep on RK4 — error message names the offending stage
// =============================================================================

namespace {

// Stateful RHS that returns the correct grid for the first `n_good`
// calls, then a wrong grid on every subsequent call. Lets us exercise
// each RK4 stage's check independently.
struct FailOnNthCall {
  Grid wrong_grid{4, 4, 4, Vec3d(0.1, 0.1, 0.1)};
  int n_good;
  mutable int calls = 0;
  Field<double> operator()(const Field<double>& psi) const {
    ++calls;
    if (calls <= n_good) return Field<double>(psi.getGrid(), 0.0);
    return Field<double>(wrong_grid, 0.0);
  }
};

void expectRk4StageInMessage(int n_good_calls,
                             const std::string& expected_stage_label) {
  Field<double> psi = makeStateField();
  try {
    integrate(psi, FailOnNthCall{Grid(4, 4, 4, Vec3d(0.1, 0.1, 0.1)), n_good_calls},
              0.01, 1, RK4{});
    FAIL() << "Expected throw at stage " << expected_stage_label;
  } catch (const std::invalid_argument& e) {
    EXPECT_NE(std::string(e.what()).find(expected_stage_label),
              std::string::npos)
        << "expected stage '" << expected_stage_label
        << "' in message: " << e.what();
  }
}

}  // namespace

TEST(IntegratorGridMismatchTest, RK4_StageLabel_K1) {
  expectRk4StageInMessage(0, "k1");
}

TEST(IntegratorGridMismatchTest, RK4_StageLabel_K2) {
  expectRk4StageInMessage(1, "k2");
}

TEST(IntegratorGridMismatchTest, RK4_StageLabel_K3) {
  expectRk4StageInMessage(2, "k3");
}

TEST(IntegratorGridMismatchTest, RK4_StageLabel_K4) {
  expectRk4StageInMessage(3, "k4");
}

// =============================================================================
// Happy-path smoke: grid-preserving RHS integrates without throwing
// =============================================================================

TEST(IntegratorGridMismatchTest, RK4_HappyPathSmoke) {
  Field<double> psi(Grid(8, 8, 8, Vec3d(0.1, 0.1, 0.1)), 1.0);
  EXPECT_NO_THROW(integrate(
      psi,
      [](const Field<double>& f) { return gridcalc::diff::laplacian(f); },
      0.0001, 3, RK4{}));
}

TEST(IntegratorGridMismatchTest, ExplicitEuler_HappyPathSmoke) {
  Field<double> psi(Grid(8, 8, 8, Vec3d(0.1, 0.1, 0.1)), 1.0);
  EXPECT_NO_THROW(integrate(
      psi,
      [](const Field<double>& f) { return gridcalc::diff::laplacian(f); },
      0.0001, 3, ExplicitEuler{}));
}
