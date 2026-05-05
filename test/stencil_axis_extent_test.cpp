// Tests for the per-axis stencil-radius precondition introduced at 0.14.1.
//
// Background: centered finite-difference stencils of radius `r` reach
// offsets [-r, r] along each active axis. On a periodic axis with extent
// N <= 2 * r, the wrap aliases offsets to the same grid sample and the
// advertised stencil silently degrades. Every diff:: and fvm:: operator
// that reads off-axis samples now rejects such grids at the operator
// entry point with std::invalid_argument.
//
// Coverage in this file:
//  - throw tests at N = 2 * radius (the largest aliasing case) for each
//    radius tier (1, 2, 3) across every affected operator;
//  - per-axis label sweep verifying the helper reports the right axis
//    (x, y, z) in its message;
//  - smallest-valid-N boundary tests (N = 2 * radius + 1) for every
//    Phase 7 / Phase 8 / Phase 13 / Phase 14 operator -- assert no-throw
//    and that the output is finite. This is the "green-side boundary"
//    coverage the FD-FFT cross-check sweep cannot give us: that sweep's
//    slope analysis is unreliable at very small N (h = O(1)), so we
//    exercise the boundary as a no-throw + finiteness gate rather than
//    by extending the slope sweep.
//  - inactive-axis carve-out verifying that an axis with per-axis
//    derivative count 0 (mixedDerivative path) accepts any positive N
//    even below the active axes' threshold.

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>

#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/diff/Biharmonic.hpp>
#include <gridcalc/diff/Divergence.hpp>
#include <gridcalc/diff/Gradient.hpp>
#include <gridcalc/diff/Laplacian.hpp>
#include <gridcalc/diff/detail/MultiIndexDerivative.hpp>
#include <gridcalc/fvm/CellLaplacian.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Mat3d;
using gridcalc::core::Vec3d;

namespace {

constexpr double kPi = 3.14159265358979323846;

Field<double> makeScalarField(int Nx, int Ny, int Nz) {
  Grid g(Nx, Ny, Nz, Vec3d(1.0, 1.0, 1.0));
  return Field<double>(g);
}

Field<Vec3d> makeVectorField(int Nx, int Ny, int Nz) {
  Grid g(Nx, Ny, Nz, Vec3d(1.0, 1.0, 1.0));
  return Field<Vec3d>(g);
}

Field<Mat3d> makeTensorField(int Nx, int Ny, int Nz) {
  Grid g(Nx, Ny, Nz, Vec3d(1.0, 1.0, 1.0));
  return Field<Mat3d>(g);
}

// Smooth single-mode initial condition on [0, 2π)^3 at extent N.
Field<double> makeSmoothScalar(int N) {
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  Grid g(N, N, N, Vec3d(h, h, h));
  return Field<double>(g, [](double x, double y, double z) {
    return std::sin(x + y + z);
  });
}

Field<Vec3d> makeSmoothVector(int N) {
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  Grid g(N, N, N, Vec3d(h, h, h));
  return Field<Vec3d>(g, [](double x, double y, double z) {
    const double s = std::sin(x + y + z);
    return Vec3d(s, s, s);
  });
}

Field<Mat3d> makeSmoothTensor(int N) {
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  Grid g(N, N, N, Vec3d(h, h, h));
  return Field<Mat3d>(g, [](double x, double y, double z) {
    Mat3d M;
    const double s = std::sin(x + y + z);
    M << s, s, s, s, s, s, s, s, s;
    return M;
  });
}

template <typename FieldT>
bool allFinite(const FieldT& f) {
  for (std::size_t n = 0; n < f.getSize(); ++n) {
    if constexpr (std::is_same_v<typename FieldT::value_type, double>) {
      if (!std::isfinite(f.data()[n])) return false;
    } else {
      // Vec3d / Mat3d: every component must be finite.
      const auto& v = f.data()[n];
      if (!v.allFinite()) return false;
    }
  }
  return true;
}

}  // namespace

// =============================================================================
// THROW TESTS — radius 1 (Order=2 Laplacian / gradient / divergence; FVM)
// =============================================================================

TEST(StencilAxisExtentTest, Radius1_LaplacianRejectsN2) {
  const auto psi = makeScalarField(2, 8, 8);
  EXPECT_THROW(gridcalc::diff::laplacian(psi), std::invalid_argument);
}

TEST(StencilAxisExtentTest, Radius1_GradientRejectsN2) {
  const auto psi = makeScalarField(2, 8, 8);
  EXPECT_THROW(gridcalc::diff::gradient(psi), std::invalid_argument);
}

TEST(StencilAxisExtentTest, Radius1_DivergenceRejectsN2) {
  const auto vfield = makeVectorField(8, 2, 8);
  EXPECT_THROW(gridcalc::diff::divergence(vfield), std::invalid_argument);
}

TEST(StencilAxisExtentTest, Radius1_TensorGradientRejectsN2) {
  const auto vfield = makeVectorField(8, 8, 2);
  EXPECT_THROW(gridcalc::diff::gradient(vfield), std::invalid_argument);
}

TEST(StencilAxisExtentTest, Radius1_TensorDivergenceRejectsN2) {
  const auto mfield = makeTensorField(2, 8, 8);
  EXPECT_THROW(gridcalc::diff::divergence(mfield), std::invalid_argument);
}

TEST(StencilAxisExtentTest, Radius1_FvmCellLaplacianRejectsN2) {
  const auto psi = makeScalarField(8, 8, 2);
  Field<double> D(psi.getGrid());
  for (std::size_t n = 0; n < D.getSize(); ++n) D.data()[n] = 1.0;
  EXPECT_THROW(gridcalc::fvm::cellLaplacian(psi, D), std::invalid_argument);
}

// =============================================================================
// THROW TESTS — radius 2 (Order=4 Laplacian / gradient / divergence;
//                         Order=2 third / fourth derivatives; biharmonic<2>)
// =============================================================================

TEST(StencilAxisExtentTest, Radius2_LaplacianOrder4RejectsN4) {
  const auto psi = makeScalarField(4, 8, 8);
  EXPECT_THROW(gridcalc::diff::laplacian<4>(psi), std::invalid_argument);
}

TEST(StencilAxisExtentTest, Radius2_GradientOrder4RejectsN4) {
  const auto psi = makeScalarField(8, 4, 8);
  EXPECT_THROW(gridcalc::diff::gradient<4>(psi), std::invalid_argument);
}

TEST(StencilAxisExtentTest, Radius2_DivergenceOrder4RejectsN4) {
  const auto vfield = makeVectorField(8, 8, 4);
  EXPECT_THROW(gridcalc::diff::divergence<4>(vfield), std::invalid_argument);
}

TEST(StencilAxisExtentTest, Radius2_BiharmonicOrder2RejectsN4) {
  const auto psi = makeScalarField(4, 8, 8);
  EXPECT_THROW(gridcalc::diff::biharmonic(psi), std::invalid_argument);
}

TEST(StencilAxisExtentTest, Radius2_MixedThirdOrder2RejectsN4) {
  // d^3/dx^3 at Order=2 -- AxisStencil<3, 2>::radius == 2.
  const auto psi = makeScalarField(4, 8, 8);
  EXPECT_THROW((gridcalc::diff::detail::mixedDerivative<3, 0, 0, 2>(psi)),
               std::invalid_argument);
}

TEST(StencilAxisExtentTest, Radius2_MixedFourthOrder2RejectsN4) {
  // d^4/dx^4 at Order=2 -- AxisStencil<4, 2>::radius == 2.
  const auto psi = makeScalarField(4, 8, 8);
  EXPECT_THROW((gridcalc::diff::detail::mixedDerivative<4, 0, 0, 2>(psi)),
               std::invalid_argument);
}

// =============================================================================
// THROW TESTS — radius 3 (Order=4 third / fourth derivatives; biharmonic<4>)
// =============================================================================

TEST(StencilAxisExtentTest, Radius3_BiharmonicOrder4RejectsN6) {
  const auto psi = makeScalarField(6, 8, 8);
  EXPECT_THROW(gridcalc::diff::biharmonic<4>(psi), std::invalid_argument);
}

TEST(StencilAxisExtentTest, Radius3_MixedThirdOrder4RejectsN6) {
  // d^3/dy^3 at Order=4 -- AxisStencil<3, 4>::radius == 3.
  const auto psi = makeScalarField(8, 6, 8);
  EXPECT_THROW((gridcalc::diff::detail::mixedDerivative<0, 3, 0, 4>(psi)),
               std::invalid_argument);
}

TEST(StencilAxisExtentTest, Radius3_MixedFourthOrder4RejectsN6) {
  // d^4/dz^4 at Order=4 -- AxisStencil<4, 4>::radius == 3.
  const auto psi = makeScalarField(8, 8, 6);
  EXPECT_THROW((gridcalc::diff::detail::mixedDerivative<0, 0, 4, 4>(psi)),
               std::invalid_argument);
}

// =============================================================================
// AXIS LABEL — each axis can independently trip the check
// =============================================================================

TEST(StencilAxisExtentTest, AxisLabel_X) {
  const auto psi = makeScalarField(4, 32, 32);
  try {
    gridcalc::diff::laplacian<4>(psi);
    FAIL() << "Expected throw";
  } catch (const std::invalid_argument& e) {
    EXPECT_NE(std::string(e.what()).find("axis x"), std::string::npos)
        << "message: " << e.what();
  }
}

TEST(StencilAxisExtentTest, AxisLabel_Y) {
  const auto psi = makeScalarField(32, 4, 32);
  try {
    gridcalc::diff::laplacian<4>(psi);
    FAIL() << "Expected throw";
  } catch (const std::invalid_argument& e) {
    EXPECT_NE(std::string(e.what()).find("axis y"), std::string::npos)
        << "message: " << e.what();
  }
}

TEST(StencilAxisExtentTest, AxisLabel_Z) {
  const auto psi = makeScalarField(32, 32, 4);
  try {
    gridcalc::diff::laplacian<4>(psi);
    FAIL() << "Expected throw";
  } catch (const std::invalid_argument& e) {
    EXPECT_NE(std::string(e.what()).find("axis z"), std::string::npos)
        << "message: " << e.what();
  }
}

// =============================================================================
// SMALLEST VALID N — green-side boundary (no throw + finite output)
// =============================================================================
//
// Each operator is exercised at N = 2 * radius + 1 (the smallest valid
// extent under the new precondition). The slope-based FD-FFT cross-check
// fixture cannot meaningfully extend its sweep to these small N values
// because h = O(1) at that resolution and the truncation-error analysis
// breaks down -- so we exercise the boundary as a no-throw + finiteness
// gate. Numerical correctness at production-relevant N is already
// verified by the existing FD-FFT sweep at {16, 32, 64, 128}.

TEST(StencilAxisExtentTest, SmallestValidN_LaplacianOrder2_NoThrow) {
  const auto psi = makeSmoothScalar(3);  // radius 1 -> N = 3
  ASSERT_NO_THROW({
    auto result = gridcalc::diff::laplacian(psi);
    EXPECT_TRUE(allFinite(result));
  });
}

TEST(StencilAxisExtentTest, SmallestValidN_LaplacianOrder4_NoThrow) {
  const auto psi = makeSmoothScalar(5);  // radius 2 -> N = 5
  ASSERT_NO_THROW({
    auto result = gridcalc::diff::laplacian<4>(psi);
    EXPECT_TRUE(allFinite(result));
  });
}

TEST(StencilAxisExtentTest, SmallestValidN_GradientOrder4_NoThrow) {
  const auto psi = makeSmoothScalar(5);
  ASSERT_NO_THROW({
    auto result = gridcalc::diff::gradient<4>(psi);
    EXPECT_TRUE(allFinite(result));
  });
}

TEST(StencilAxisExtentTest, SmallestValidN_DivergenceOrder4_NoThrow) {
  const auto vfield = makeSmoothVector(5);
  ASSERT_NO_THROW({
    auto result = gridcalc::diff::divergence<4>(vfield);
    EXPECT_TRUE(allFinite(result));
  });
}

TEST(StencilAxisExtentTest, SmallestValidN_TensorGradientOrder4_NoThrow) {
  const auto vfield = makeSmoothVector(5);
  ASSERT_NO_THROW({
    auto result = gridcalc::diff::gradient<4>(vfield);
    EXPECT_TRUE(allFinite(result));
  });
}

TEST(StencilAxisExtentTest, SmallestValidN_TensorDivergenceOrder4_NoThrow) {
  const auto mfield = makeSmoothTensor(5);
  ASSERT_NO_THROW({
    auto result = gridcalc::diff::divergence<4>(mfield);
    EXPECT_TRUE(allFinite(result));
  });
}

TEST(StencilAxisExtentTest, SmallestValidN_BiharmonicOrder2_NoThrow) {
  const auto psi = makeSmoothScalar(5);  // max(S4::radius=2, S2::radius=1) = 2
  ASSERT_NO_THROW({
    auto result = gridcalc::diff::biharmonic(psi);
    EXPECT_TRUE(allFinite(result));
  });
}

TEST(StencilAxisExtentTest, SmallestValidN_BiharmonicOrder4_NoThrow) {
  const auto psi = makeSmoothScalar(7);  // max(S4::radius=3, S2::radius=2) = 3
  ASSERT_NO_THROW({
    auto result = gridcalc::diff::biharmonic<4>(psi);
    EXPECT_TRUE(allFinite(result));
  });
}

TEST(StencilAxisExtentTest, SmallestValidN_MixedThirdOrder2_NoThrow) {
  const auto psi = makeSmoothScalar(5);  // radius 2 -> N = 5
  try {
    auto result = gridcalc::diff::detail::mixedDerivative<3, 0, 0, 2>(psi);
    EXPECT_TRUE(allFinite(result));
  } catch (...) {
    FAIL() << "Unexpected throw at smallest valid N";
  }
}

TEST(StencilAxisExtentTest, SmallestValidN_MixedFourthOrder4_NoThrow) {
  const auto psi = makeSmoothScalar(7);  // radius 3 -> N = 7
  try {
    auto result = gridcalc::diff::detail::mixedDerivative<4, 0, 0, 4>(psi);
    EXPECT_TRUE(allFinite(result));
  } catch (...) {
    FAIL() << "Unexpected throw at smallest valid N";
  }
}

TEST(StencilAxisExtentTest, SmallestValidN_FvmCellLaplacian_NoThrow) {
  const auto psi = makeSmoothScalar(3);
  Field<double> D(psi.getGrid());
  for (std::size_t n = 0; n < D.getSize(); ++n) D.data()[n] = 1.0;
  ASSERT_NO_THROW({
    auto result = gridcalc::fvm::cellLaplacian(psi, D);
    EXPECT_TRUE(allFinite(result));
  });
}

// =============================================================================
// DEGENERATE AXIS (N = 1) — every offset wraps to index 0; the centered
// stencil weights sum to zero by construction, so the aliased result is
// exactly 0 (the correct value for a single-cell axis with no spatial
// variation). Accepted for every operator regardless of stencil radius.
// =============================================================================

TEST(StencilAxisExtentTest, DegenerateAxis_LaplacianOrder4_AcceptsN1) {
  const auto psi = makeScalarField(8, 8, 1);
  Field<double> result(psi.getGrid());
  ASSERT_NO_THROW(result = gridcalc::diff::laplacian<4>(psi));
  EXPECT_TRUE(allFinite(result));
}

TEST(StencilAxisExtentTest, DegenerateAxis_GradientOrder4_AcceptsN1) {
  const auto psi = makeScalarField(8, 1, 8);
  Field<Vec3d> result(psi.getGrid());
  ASSERT_NO_THROW(result = gridcalc::diff::gradient<4>(psi));
  EXPECT_TRUE(allFinite(result));
}

TEST(StencilAxisExtentTest, DegenerateAxis_BiharmonicOrder4_AcceptsN1) {
  // N=1 along z; biharmonic still reads radius-3 along z but sums to 0.
  const auto psi = makeScalarField(8, 8, 1);
  Field<double> result(psi.getGrid());
  ASSERT_NO_THROW(result = gridcalc::diff::biharmonic<4>(psi));
  EXPECT_TRUE(allFinite(result));
}

TEST(StencilAxisExtentTest, DegenerateAxis_FvmCellLaplacian_AcceptsN1) {
  // 1D-flavored FVM use case (Nx active, Ny=Nz=1 degenerate). Documents
  // the existing FvmCellLaplacianTest.HarmonicMeanFaceCoefficient
  // pattern, which would otherwise be broken by a strict
  // N >= 2*radius+1 check.
  const auto psi = makeScalarField(4, 1, 1);
  Field<double> D(psi.getGrid());
  for (std::size_t n = 0; n < D.getSize(); ++n) D.data()[n] = 1.0;
  Field<double> result(psi.getGrid());
  ASSERT_NO_THROW(result = gridcalc::fvm::cellLaplacian(psi, D));
  EXPECT_TRUE(allFinite(result));
}

// =============================================================================
// INACTIVE AXIS — N_a = 0 carve-out: tiny extent on an unused axis is fine
// =============================================================================

TEST(StencilAxisExtentTest, InactiveAxisAcceptsTinyExtent) {
  // d^3/dx^3 at Order=2: x has radius 2, y and z have radius 0
  // (AxisStencil<0, *>) so N_y = N_z = 1 is accepted.
  const auto psi = makeScalarField(8, 1, 1);
  try {
    auto result = gridcalc::diff::detail::mixedDerivative<3, 0, 0, 2>(psi);
    EXPECT_TRUE(allFinite(result));
  } catch (...) {
    FAIL() << "Unexpected throw on inactive-axis path";
  }
}
