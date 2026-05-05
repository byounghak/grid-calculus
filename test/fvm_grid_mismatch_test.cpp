// Tests for the 0.14.5 D-grid + D > 0 preconditions on the FVM
// heterogeneous-D path, plus the CFL-tightening from D_max to
// D_max_face.
//
// Coverage:
//  - fvm::cellLaplacian rejects on D-grid mismatch (three flavors:
//    different total cell count, same total / different per-axis
//    extents, same shape / different cell size);
//  - fvm::cellLaplacian rejects on D <= 0 at any cell, with a flat-
//    index message;
//  - solve::diffuse(psi, D, ...) rejects on D-grid mismatch with a
//    "solve::diffuse: D" context (catches the bug before the
//    integrator starts);
//  - solve::diffuse(psi, D, ...) rejects on D <= 0 (now via
//    cellLaplacian on the first integrator stage);
//  - happy-path smoke: a correctly-grid-D field integrates one step;
//  - CFL tightening on a high-contrast D field: a dt that exceeded the
//    old D_max bound is accepted under the new D_max_face bound;
//  - uniform-D regression: D_max_face = D_max for uniform D, so the
//    new bound is bit-identical to the old.

#include <gtest/gtest.h>

#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/fvm/CellLaplacian.hpp>
#include <gridcalc/solve/Diffusion.hpp>
#include <gridcalc/solve/ExplicitEuler.hpp>
#include <gridcalc/solve/RK4.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
using gridcalc::solve::ExplicitEuler;
using gridcalc::solve::RK4;
namespace gs = gridcalc::solve;

namespace {

Field<double> makeUnitFieldOnGrid(const Grid& g, double v) {
  Field<double> f(g);
  for (std::size_t n = 0; n < f.getSize(); ++n) f.data()[n] = v;
  return f;
}

Grid makeStandardGrid() {
  return Grid(8, 8, 8, Vec3d(0.1, 0.1, 0.1));
}

}  // namespace

// =============================================================================
// fvm::cellLaplacian — D-grid mismatch (three flavors)
// =============================================================================

TEST(FvmGridMismatchTest, CellLaplacian_RejectsDifferentTotalCellCount) {
  const Grid g_psi = makeStandardGrid();         // 8x8x8 = 512
  const Grid g_D(4, 4, 4, Vec3d(0.1, 0.1, 0.1));  // 64
  const auto psi = makeUnitFieldOnGrid(g_psi, 0.0);
  const auto D = makeUnitFieldOnGrid(g_D, 1.0);
  EXPECT_THROW(gridcalc::fvm::cellLaplacian(psi, D), std::invalid_argument);
}

TEST(FvmGridMismatchTest, CellLaplacian_RejectsSameSizeDifferentShape) {
  const Grid g_psi = makeStandardGrid();           // 8x8x8 = 512
  const Grid g_D(4, 16, 8, Vec3d(0.1, 0.1, 0.1));  // 4*16*8 = 512
  const auto psi = makeUnitFieldOnGrid(g_psi, 0.0);
  const auto D = makeUnitFieldOnGrid(g_D, 1.0);
  EXPECT_THROW(gridcalc::fvm::cellLaplacian(psi, D), std::invalid_argument);
}

TEST(FvmGridMismatchTest, CellLaplacian_RejectsSameShapeDifferentCellSize) {
  const Grid g_psi = makeStandardGrid();           // h = 0.1
  const Grid g_D(8, 8, 8, Vec3d(0.2, 0.2, 0.2));   // h = 0.2
  const auto psi = makeUnitFieldOnGrid(g_psi, 0.0);
  const auto D = makeUnitFieldOnGrid(g_D, 1.0);
  EXPECT_THROW(gridcalc::fvm::cellLaplacian(psi, D), std::invalid_argument);
}

// =============================================================================
// fvm::cellLaplacian — D > 0 violation
// =============================================================================

TEST(FvmGridMismatchTest, CellLaplacian_RejectsZeroD) {
  const Grid g = makeStandardGrid();
  const auto psi = makeUnitFieldOnGrid(g, 0.0);
  Field<double> D = makeUnitFieldOnGrid(g, 1.0);
  D(2, 3, 4) = 0.0;  // single zero cell
  try {
    gridcalc::fvm::cellLaplacian(psi, D);
    FAIL() << "Expected throw";
  } catch (const std::invalid_argument& e) {
    const std::string msg(e.what());
    EXPECT_NE(msg.find("flat index"), std::string::npos)
        << "message should reference flat index: " << msg;
  }
}

TEST(FvmGridMismatchTest, CellLaplacian_RejectsNegativeD) {
  const Grid g = makeStandardGrid();
  const auto psi = makeUnitFieldOnGrid(g, 0.0);
  Field<double> D = makeUnitFieldOnGrid(g, 1.0);
  D(0, 0, 0) = -0.5;
  EXPECT_THROW(gridcalc::fvm::cellLaplacian(psi, D), std::invalid_argument);
}

// =============================================================================
// solve::diffuse — D-grid mismatch (three flavors), with clearer context
// =============================================================================

TEST(FvmGridMismatchTest, Diffuse_RejectsDifferentTotalCellCount) {
  const Grid g_psi = makeStandardGrid();
  const Grid g_D(4, 4, 4, Vec3d(0.1, 0.1, 0.1));
  Field<double> psi = makeUnitFieldOnGrid(g_psi, 1.0);
  const auto D = makeUnitFieldOnGrid(g_D, 1.0);
  EXPECT_THROW(gs::diffuse(psi, D, 0.001, 1, ExplicitEuler{}),
               std::invalid_argument);
}

TEST(FvmGridMismatchTest, Diffuse_RejectsSameSizeDifferentShape) {
  const Grid g_psi = makeStandardGrid();
  const Grid g_D(4, 16, 8, Vec3d(0.1, 0.1, 0.1));
  Field<double> psi = makeUnitFieldOnGrid(g_psi, 1.0);
  const auto D = makeUnitFieldOnGrid(g_D, 1.0);
  EXPECT_THROW(gs::diffuse(psi, D, 0.001, 1, ExplicitEuler{}),
               std::invalid_argument);
}

TEST(FvmGridMismatchTest, Diffuse_RejectsSameShapeDifferentCellSize) {
  const Grid g_psi = makeStandardGrid();
  const Grid g_D(8, 8, 8, Vec3d(0.2, 0.2, 0.2));
  Field<double> psi = makeUnitFieldOnGrid(g_psi, 1.0);
  const auto D = makeUnitFieldOnGrid(g_D, 1.0);
  EXPECT_THROW(gs::diffuse(psi, D, 0.001, 1, ExplicitEuler{}),
               std::invalid_argument);
}

TEST(FvmGridMismatchTest, Diffuse_RejectsDGridMismatch_NamesContext) {
  const Grid g_psi = makeStandardGrid();
  const Grid g_D(4, 4, 4, Vec3d(0.1, 0.1, 0.1));
  Field<double> psi = makeUnitFieldOnGrid(g_psi, 1.0);
  const auto D = makeUnitFieldOnGrid(g_D, 1.0);
  try {
    gs::diffuse(psi, D, 0.001, 1, ExplicitEuler{});
    FAIL() << "Expected throw";
  } catch (const std::invalid_argument& e) {
    const std::string msg(e.what());
    EXPECT_NE(msg.find("solve::diffuse"), std::string::npos)
        << "message should name the solve::diffuse entry point: " << msg;
  }
}

// =============================================================================
// solve::diffuse — D <= 0 propagated from cellLaplacian
// =============================================================================

TEST(FvmGridMismatchTest, Diffuse_RejectsZeroDViaCellLaplacian) {
  const Grid g = makeStandardGrid();
  Field<double> psi = makeUnitFieldOnGrid(g, 1.0);
  Field<double> D = makeUnitFieldOnGrid(g, 1.0);
  D(0, 0, 0) = 0.0;
  EXPECT_THROW(gs::diffuse(psi, D, 0.001, 1, ExplicitEuler{}),
               std::invalid_argument);
}

// =============================================================================
// Happy-path smoke: a correctly-grid-D field integrates without throwing
// =============================================================================

TEST(FvmGridMismatchTest, CellLaplacian_HappyPathSmoke) {
  const Grid g = makeStandardGrid();
  Field<double> psi(g, [](double x, double, double) { return std::sin(x); });
  Field<double> D = makeUnitFieldOnGrid(g, 1.0);
  EXPECT_NO_THROW(gridcalc::fvm::cellLaplacian(psi, D));
}

TEST(FvmGridMismatchTest, Diffuse_HappyPathSmoke) {
  const Grid g = makeStandardGrid();
  Field<double> psi(g, [](double x, double, double) { return std::sin(x); });
  Field<double> D = makeUnitFieldOnGrid(g, 1.0);
  // dt small enough for explicit Euler with D = 1 on h = 0.1 grid:
  //   D * dt * 3/h^2 <= 0.5  =>  dt <= 0.5 * 0.01 / 3 = 1.67e-3.
  EXPECT_NO_THROW(gs::diffuse(psi, D, 1.0e-4, 5, ExplicitEuler{}));
}

// =============================================================================
// CFL tightening — high-contrast D admits dt that would have thrown
// =============================================================================
//
// Build D = 1 everywhere except one hot cell at D = 100.
//   harmonicMean(1, 100) = 200/101 ≈ 1.98
//   D_max_face ≈ 1.98 (at the 6 faces of the hot cell)
//   D_max     = 100 (at the hot cell itself)
//
// Pick dt such that:
//   dt * 100   * 3/h^2 > 0.5  (would throw under the old bound)
//   dt *   1.98 * 3/h^2 < 0.5  (passes under the new bound)
//
// h = 0.1, sum_a (1/h^2) = 300 (cubic).
//   dt > 0.5 / (100 * 300) = 1.667e-5
//   dt < 0.5 / (1.98 * 300) = 8.418e-4
//
// dt = 1e-4 is comfortably between these bounds.

TEST(FvmGridMismatchTest, Diffuse_HighContrastD_AdmitsLooserDt) {
  const Grid g = makeStandardGrid();
  Field<double> psi = makeUnitFieldOnGrid(g, 1.0);
  Field<double> D = makeUnitFieldOnGrid(g, 1.0);
  D(4, 4, 4) = 100.0;  // single hot cell
  const double dt = 1.0e-4;
  // Old bound: dt * 100 * 300 = 3.0  -- way over 0.5, would have thrown.
  // New bound: dt * 1.98 * 300 ~ 0.0594  -- well under 0.5, passes.
  EXPECT_NO_THROW(gs::diffuse(psi, D, dt, 1, ExplicitEuler{}));
}

TEST(FvmGridMismatchTest, Diffuse_UniformD_BitIdenticalCFLBehavior) {
  // For uniform D, harmonicMean(D, D) = D, so D_max_face = D_max.
  // The CFL boundary is at the same dt as before. Pick a dt right at
  // the edge: D * dt * 3/h^2 = 0.5 - eps should pass; +eps should
  // throw. This regression test confirms no behavioral change on
  // uniform D.
  const Grid g = makeStandardGrid();
  Field<double> psi = makeUnitFieldOnGrid(g, 1.0);
  Field<double> D = makeUnitFieldOnGrid(g, 1.0);
  // sum(1/h^2) = 300; D = 1; CFL limit = 0.5; boundary at dt = 0.5/300
  //   = 1.667e-3. Just under: 1.6e-3 should pass.
  EXPECT_NO_THROW(gs::diffuse(psi, D, 1.6e-3, 1, ExplicitEuler{}));
  // Just over: 1.8e-3 should throw on the CFL gate.
  Field<double> psi2 = makeUnitFieldOnGrid(g, 1.0);
  EXPECT_THROW(gs::diffuse(psi2, D, 1.8e-3, 1, ExplicitEuler{}),
               std::invalid_argument);
}
