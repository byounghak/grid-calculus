#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/diff/Biharmonic.hpp>
#include <gridcalc/diff/FourthOrder.hpp>
#include <gridcalc/diff/MixedPartial.hpp>
#include <gridcalc/diff/ThirdOrder.hpp>
#include <gridcalc/diff/detail/MultiIndexDerivative.hpp>
#include <gridcalc/stencil/FourthDerivative.hpp>
#include <gridcalc/stencil/ThirdDerivative.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
using gridcalc::stencil::FourthDerivative;
using gridcalc::stencil::ThirdDerivative;

namespace gd = gridcalc::diff;

namespace {

constexpr double kPi = 3.14159265358979323846;

double maxAbs(const Field<double>& f) {
  double m = 0.0;
  for (double v : f) {
    const double a = std::abs(v);
    if (a > m) m = a;
  }
  return m;
}

double maxAbsDiff(const Field<double>& a, const Field<double>& b) {
  double m = 0.0;
  for (std::size_t n = 0; n < a.getSize(); ++n) {
    const double d = std::abs(a.data()[n] - b.data()[n]);
    if (d > m) m = d;
  }
  return m;
}

double logLogSlope(const std::vector<int>& Ns, const std::vector<double>& errors,
                   double L) {
  const std::size_t n = Ns.size();
  double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double h = L / static_cast<double>(Ns[i]);
    const double x = std::log(h);
    const double y = std::log(errors[i]);
    sum_x += x;
    sum_y += y;
    sum_xy += x * y;
    sum_xx += x * x;
  }
  return (static_cast<double>(n) * sum_xy - sum_x * sum_y) /
         (static_cast<double>(n) * sum_xx - sum_x * sum_x);
}

// d^n / dx^n sin(x). 4-cycle: sin, cos, -sin, -cos.
double sinDerivative(int n, double x) {
  switch (n % 4) {
    case 0: return std::sin(x);
    case 1: return std::cos(x);
    case 2: return -std::sin(x);
    case 3: return -std::cos(x);
  }
  return 0.0;
}

// Reference: ∂^(Nx+Ny+Nz)/∂x^Nx ∂y^Ny ∂z^Nz of sin(x)sin(y)sin(z).
template <int Nx, int Ny, int Nz>
double analyticalPartial(double x, double y, double z) {
  return sinDerivative(Nx, x) * sinDerivative(Ny, y) * sinDerivative(Nz, z);
}

// Run the convergence sweep on N ∈ {16, 32, 64} for the given multi-index
// using the provided callable `op(field) -> Field<double>`. Returns the
// log-log slope of the relative max-norm error vs h.
template <int Nx, int Ny, int Nz, typename Op>
double sweepSlope(Op&& op) {
  const double L = 2.0 * kPi;
  const std::vector<int> Ns = {16, 32, 64};
  std::vector<double> errors;
  for (int N : Ns) {
    const double h = L / static_cast<double>(N);
    Grid grid(N, N, N, Vec3d(h, h, h));
    Field<double> psi(grid, [](double x, double y, double z) {
      return std::sin(x) * std::sin(y) * std::sin(z);
    });
    Field<double> expected(grid, [](double x, double y, double z) {
      return analyticalPartial<Nx, Ny, Nz>(x, y, z);
    });
    Field<double> result = op(psi);
    const double ref = std::max(maxAbs(expected), 1e-300);
    errors.push_back(maxAbsDiff(result, expected) / ref);
  }
  return logLogSlope(Ns, errors, L);
}

}  // namespace

// =============================================================================
// Weight-table audits
// =============================================================================

TEST(Phase8WeightsTest, ThirdDerivativeOrder2) {
  EXPECT_EQ(ThirdDerivative<2>::radius, 2);
  EXPECT_DOUBLE_EQ(ThirdDerivative<2>::weights[0], -1.0 / 2.0);
  EXPECT_DOUBLE_EQ(ThirdDerivative<2>::weights[1], 1.0);
  EXPECT_DOUBLE_EQ(ThirdDerivative<2>::weights[2], 0.0);
  EXPECT_DOUBLE_EQ(ThirdDerivative<2>::weights[3], -1.0);
  EXPECT_DOUBLE_EQ(ThirdDerivative<2>::weights[4], 1.0 / 2.0);
}

TEST(Phase8WeightsTest, ThirdDerivativeOrder4) {
  EXPECT_EQ(ThirdDerivative<4>::radius, 3);
  EXPECT_DOUBLE_EQ(ThirdDerivative<4>::weights[0], 1.0 / 8.0);
  EXPECT_DOUBLE_EQ(ThirdDerivative<4>::weights[1], -1.0);
  EXPECT_DOUBLE_EQ(ThirdDerivative<4>::weights[2], 13.0 / 8.0);
  EXPECT_DOUBLE_EQ(ThirdDerivative<4>::weights[3], 0.0);
  EXPECT_DOUBLE_EQ(ThirdDerivative<4>::weights[4], -13.0 / 8.0);
  EXPECT_DOUBLE_EQ(ThirdDerivative<4>::weights[5], 1.0);
  EXPECT_DOUBLE_EQ(ThirdDerivative<4>::weights[6], -1.0 / 8.0);
}

TEST(Phase8WeightsTest, FourthDerivativeOrder2) {
  EXPECT_EQ(FourthDerivative<2>::radius, 2);
  EXPECT_DOUBLE_EQ(FourthDerivative<2>::weights[0], 1.0);
  EXPECT_DOUBLE_EQ(FourthDerivative<2>::weights[1], -4.0);
  EXPECT_DOUBLE_EQ(FourthDerivative<2>::weights[2], 6.0);
  EXPECT_DOUBLE_EQ(FourthDerivative<2>::weights[3], -4.0);
  EXPECT_DOUBLE_EQ(FourthDerivative<2>::weights[4], 1.0);
}

TEST(Phase8WeightsTest, FourthDerivativeOrder4) {
  EXPECT_EQ(FourthDerivative<4>::radius, 3);
  EXPECT_DOUBLE_EQ(FourthDerivative<4>::weights[0], -1.0 / 6.0);
  EXPECT_DOUBLE_EQ(FourthDerivative<4>::weights[1], 2.0);
  EXPECT_DOUBLE_EQ(FourthDerivative<4>::weights[2], -13.0 / 2.0);
  EXPECT_DOUBLE_EQ(FourthDerivative<4>::weights[3], 28.0 / 3.0);
  EXPECT_DOUBLE_EQ(FourthDerivative<4>::weights[4], -13.0 / 2.0);
  EXPECT_DOUBLE_EQ(FourthDerivative<4>::weights[5], 2.0);
  EXPECT_DOUBLE_EQ(FourthDerivative<4>::weights[6], -1.0 / 6.0);
}

// =============================================================================
// Convergence sweeps for every named partial-derivative function
//
// Each TEST runs the N ∈ {16, 32, 64} sweep on the trigonometric manufactured
// solution ψ = sin(x)sin(y)sin(z) over the periodic [0, 2π)^3 box and asserts
// the log-log slope of relative max-norm error vs h falls within
// [Order − 0.5, Order + 0.5]. The named function is exercised explicitly so
// the API binding is part of the test.
// =============================================================================

#define EXPECT_PARTIAL_SLOPE_RANGE(NAME, NX, NY, NZ, ORDER)                  \
  TEST(Phase8ConvergenceTest, NAME##_Order##ORDER##_Sweep) {                  \
    const double slope = sweepSlope<NX, NY, NZ>(                              \
        [](const Field<double>& f) { return gd::NAME<ORDER>(f); });           \
    EXPECT_GT(slope, ORDER - 0.5) << #NAME "<" #ORDER "> slope=" << slope;    \
    EXPECT_LT(slope, ORDER + 0.5) << #NAME "<" #ORDER "> slope=" << slope;    \
  }

// Rank-2 partials.
EXPECT_PARTIAL_SLOPE_RANGE(dxx, 2, 0, 0, 2)
EXPECT_PARTIAL_SLOPE_RANGE(dxx, 2, 0, 0, 4)
EXPECT_PARTIAL_SLOPE_RANGE(dyy, 0, 2, 0, 2)
EXPECT_PARTIAL_SLOPE_RANGE(dyy, 0, 2, 0, 4)
EXPECT_PARTIAL_SLOPE_RANGE(dzz, 0, 0, 2, 2)
EXPECT_PARTIAL_SLOPE_RANGE(dzz, 0, 0, 2, 4)
EXPECT_PARTIAL_SLOPE_RANGE(dxy, 1, 1, 0, 2)
EXPECT_PARTIAL_SLOPE_RANGE(dxy, 1, 1, 0, 4)
EXPECT_PARTIAL_SLOPE_RANGE(dxz, 1, 0, 1, 2)
EXPECT_PARTIAL_SLOPE_RANGE(dxz, 1, 0, 1, 4)
EXPECT_PARTIAL_SLOPE_RANGE(dyz, 0, 1, 1, 2)
EXPECT_PARTIAL_SLOPE_RANGE(dyz, 0, 1, 1, 4)

// Rank-3 partials.
EXPECT_PARTIAL_SLOPE_RANGE(d3dx3, 3, 0, 0, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d3dx3, 3, 0, 0, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d3dy3, 0, 3, 0, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d3dy3, 0, 3, 0, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d3dz3, 0, 0, 3, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d3dz3, 0, 0, 3, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d3dx2dy, 2, 1, 0, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d3dx2dy, 2, 1, 0, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d3dxdy2, 1, 2, 0, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d3dxdy2, 1, 2, 0, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d3dx2dz, 2, 0, 1, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d3dx2dz, 2, 0, 1, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d3dxdz2, 1, 0, 2, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d3dxdz2, 1, 0, 2, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d3dy2dz, 0, 2, 1, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d3dy2dz, 0, 2, 1, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d3dydz2, 0, 1, 2, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d3dydz2, 0, 1, 2, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d3dxdydz, 1, 1, 1, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d3dxdydz, 1, 1, 1, 4)

// Rank-4 partials.
EXPECT_PARTIAL_SLOPE_RANGE(d4dx4, 4, 0, 0, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d4dx4, 4, 0, 0, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d4dy4, 0, 4, 0, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d4dy4, 0, 4, 0, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d4dz4, 0, 0, 4, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d4dz4, 0, 0, 4, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d4dx3dy, 3, 1, 0, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d4dx3dy, 3, 1, 0, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d4dxdy3, 1, 3, 0, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d4dxdy3, 1, 3, 0, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d4dx3dz, 3, 0, 1, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d4dx3dz, 3, 0, 1, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d4dxdz3, 1, 0, 3, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d4dxdz3, 1, 0, 3, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d4dy3dz, 0, 3, 1, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d4dy3dz, 0, 3, 1, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d4dydz3, 0, 1, 3, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d4dydz3, 0, 1, 3, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d4dx2dy2, 2, 2, 0, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d4dx2dy2, 2, 2, 0, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d4dx2dz2, 2, 0, 2, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d4dx2dz2, 2, 0, 2, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d4dy2dz2, 0, 2, 2, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d4dy2dz2, 0, 2, 2, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d4dx2dydz, 2, 1, 1, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d4dx2dydz, 2, 1, 1, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d4dxdy2dz, 1, 2, 1, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d4dxdy2dz, 1, 2, 1, 4)
EXPECT_PARTIAL_SLOPE_RANGE(d4dxdydz2, 1, 1, 2, 2)
EXPECT_PARTIAL_SLOPE_RANGE(d4dxdydz2, 1, 1, 2, 4)

#undef EXPECT_PARTIAL_SLOPE_RANGE

// =============================================================================
// Biharmonic: convergence sweep on ψ = sin(x)sin(y)sin(z) (∇⁴ψ = 9·ψ).
// =============================================================================

namespace {

template <int Order, typename Op>
double biharmonicSweepSlope(Op&& op, double k1, double k2, double k3) {
  const double L = 2.0 * kPi;
  const std::vector<int> Ns = {16, 32, 64};
  const double k_sq = k1 * k1 + k2 * k2 + k3 * k3;
  const double k_fourth = k_sq * k_sq;
  std::vector<double> errors;
  for (int N : Ns) {
    const double h = L / static_cast<double>(N);
    Grid grid(N, N, N, Vec3d(h, h, h));
    Field<double> psi(grid, [k1, k2, k3](double x, double y, double z) {
      return std::sin(k1 * x + k2 * y + k3 * z);
    });
    Field<double> expected(
        grid, [k1, k2, k3, k_fourth](double x, double y, double z) {
          return k_fourth * std::sin(k1 * x + k2 * y + k3 * z);
        });
    Field<double> result = op(psi);
    const double ref = std::max(maxAbs(expected), 1e-300);
    errors.push_back(maxAbsDiff(result, expected) / ref);
  }
  return logLogSlope(Ns, errors, L);
}

}  // namespace

TEST(Phase8BiharmonicTest, Order2Sweep) {
  const double slope = biharmonicSweepSlope<2>(
      [](const Field<double>& f) { return gd::biharmonic<2>(f); }, 1.0, 1.0, 1.0);
  EXPECT_GT(slope, 1.5) << "biharmonic<2> slope=" << slope;
  EXPECT_LT(slope, 2.5) << "biharmonic<2> slope=" << slope;
}

TEST(Phase8BiharmonicTest, Order4Sweep) {
  const double slope = biharmonicSweepSlope<4>(
      [](const Field<double>& f) { return gd::biharmonic<4>(f); }, 1.0, 1.0, 1.0);
  EXPECT_GT(slope, 3.5) << "biharmonic<4> slope=" << slope;
  EXPECT_LT(slope, 4.5) << "biharmonic<4> slope=" << slope;
}

// =============================================================================
// Roadmap acceptance: biharmonic of sin(k·x) recovers |k|⁴.
// At Order 4 on N = 64 the relative max-norm error against |k|⁴·ψ should sit
// well below 5% for k = (1, 1, 1) (|k|² = 3, |k|⁴ = 9) and k = (1, 2, 1)
// (|k|² = 6, |k|⁴ = 36). Order-2 tolerance is looser to match its truncation.
// =============================================================================

namespace {

double biharmonicRelativeError(int N, double k1, double k2, double k3,
                               int order) {
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [k1, k2, k3](double x, double y, double z) {
    return std::sin(k1 * x + k2 * y + k3 * z);
  });
  const double k_sq = k1 * k1 + k2 * k2 + k3 * k3;
  const double k_fourth = k_sq * k_sq;
  Field<double> expected(grid,
                         [k1, k2, k3, k_fourth](double x, double y, double z) {
                           return k_fourth *
                                  std::sin(k1 * x + k2 * y + k3 * z);
                         });
  Field<double> result =
      (order == 2) ? gd::biharmonic<2>(psi) : gd::biharmonic<4>(psi);
  return maxAbsDiff(result, expected) / std::max(maxAbs(expected), 1e-300);
}

}  // namespace

TEST(Phase8BiharmonicTest, RecoversKToTheFourthOrder2_KTriple_1_1_1) {
  // |k|² = 3; expected truncation at N=64, order 2: O(h² · |k|² · |k|⁴) /
  // |k|⁴ ~ h²·|k|² ≈ 0.029. Allow a factor of 3 safety margin.
  EXPECT_LT(biharmonicRelativeError(64, 1.0, 1.0, 1.0, 2), 0.10);
}

TEST(Phase8BiharmonicTest, RecoversKToTheFourthOrder2_KTriple_1_2_1) {
  // |k|² = 6; expected truncation at N=64, order 2 ~ h²·|k|² ≈ 0.058.
  EXPECT_LT(biharmonicRelativeError(64, 1.0, 2.0, 1.0, 2), 0.20);
}

TEST(Phase8BiharmonicTest, RecoversKToTheFourthOrder4_KTriple_1_1_1) {
  // Order 4 on N=64: relative error ~ h⁴·|k|² ≈ 2.7e-4.
  EXPECT_LT(biharmonicRelativeError(64, 1.0, 1.0, 1.0, 4), 1.0e-3);
}

TEST(Phase8BiharmonicTest, RecoversKToTheFourthOrder4_KTriple_1_2_1) {
  // Order 4 on N=64, |k|² = 6.
  EXPECT_LT(biharmonicRelativeError(64, 1.0, 2.0, 1.0, 4), 5.0e-3);
}

// =============================================================================
// d4 / biharmonic alias parity.
// =============================================================================

TEST(Phase8BiharmonicTest, D4AliasMatchesBiharmonicOrder2) {
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [](double x, double y, double z) {
    return std::sin(x) * std::sin(y) * std::sin(z);
  });
  Field<double> via_alias = gd::d4<2>(psi);
  Field<double> via_canonical = gd::biharmonic<2>(psi);
  for (std::size_t i = 0; i < psi.getSize(); ++i) {
    EXPECT_DOUBLE_EQ(via_alias.data()[i], via_canonical.data()[i]);
  }
}

TEST(Phase8BiharmonicTest, D4AliasMatchesBiharmonicOrder4) {
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [](double x, double y, double z) {
    return std::sin(x) * std::sin(y) * std::sin(z);
  });
  Field<double> via_alias = gd::d4<4>(psi);
  Field<double> via_canonical = gd::biharmonic<4>(psi);
  for (std::size_t i = 0; i < psi.getSize(); ++i) {
    EXPECT_DOUBLE_EQ(via_alias.data()[i], via_canonical.data()[i]);
  }
}
