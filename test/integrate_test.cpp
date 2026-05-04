#include <gtest/gtest.h>

#include <cmath>
#include <cstring>

#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/func/Integrate.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
using gridcalc::func::integrate;
using gridcalc::func::Kahan;
using gridcalc::func::Pairwise;

namespace {

constexpr double kPi = 3.14159265358979323846;

bool bitEqual(double a, double b) {
  // Bit-exact comparison via memcpy to avoid signed-NaN edge cases.
  std::uint64_t ua, ub;
  std::memcpy(&ua, &a, sizeof(double));
  std::memcpy(&ub, &b, sizeof(double));
  return ua == ub;
}

}  // namespace

TEST(IntegrateTest, UnitFieldRecoversDomainVolume) {
  // Anisotropic, non-cubic grid to catch any axis-mixup in the volume element.
  const int Nx = 8;
  const int Ny = 16;
  const int Nz = 4;
  const double hx = 0.25;
  const double hy = 0.125;
  const double hz = 0.5;
  const double expected_volume = (Nx * hx) * (Ny * hy) * (Nz * hz);

  Grid grid(Nx, Ny, Nz, Vec3d(hx, hy, hz));
  Field<double> ones(grid, 1.0);

  EXPECT_DOUBLE_EQ(integrate(ones), expected_volume);
  EXPECT_DOUBLE_EQ(integrate(ones, Pairwise{}), expected_volume);
  EXPECT_DOUBLE_EQ(integrate(ones, Kahan{}), expected_volume);
}

TEST(IntegrateTest, SinOverPeriodIsZero) {
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  const double kx = 1.0;

  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [&](double x, double, double) { return std::sin(kx * x); });

  // Periodic trapezoidal rule integrates a single sinusoid exactly to 0
  // (up to rounding); allow a tight tolerance.
  EXPECT_NEAR(integrate(psi), 0.0, 1e-12);
  EXPECT_NEAR(integrate(psi, Kahan{}), 0.0, 1e-12);
}

TEST(IntegrateTest, ManufacturedSolutionPeriodicProduct) {
  // f(x,y,z) = cos(x) sin(2y) sin(3z) on [0, 2pi]^3.
  // Each axis integrates to zero independently, so the integral is exactly 0
  // (up to rounding) under the periodic midpoint rule.
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);

  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> f(grid, [&](double x, double y, double z) {
    return std::cos(x) * std::sin(2.0 * y) * std::sin(3.0 * z);
  });

  EXPECT_NEAR(integrate(f), 0.0, 1e-12);
  EXPECT_NEAR(integrate(f, Kahan{}), 0.0, 1e-12);
}

TEST(IntegrateTest, PairwiseAndKahanAgree) {
  // Integrand chosen so the analytical integral over [0, 2pi]^3 is non-zero
  // (the constant term contributes (2*pi)^3; the trig term integrates to 0):
  //   f(x,y,z) = 1 + sin(x)*cos(y).
  // Discrete integral therefore sits near (2*pi)^3 ~= 248, giving a stable
  // scale for the relative-error comparison.
  const int N = 64;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);

  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> f(grid, [&](double x, double y, double) {
    return 1.0 + std::sin(x) * std::cos(y);
  });

  const double pairwise = integrate(f, Pairwise{});
  const double kahan = integrate(f, Kahan{});
  const double scale = std::max(std::abs(pairwise), std::abs(kahan));
  EXPECT_LT(std::abs(pairwise - kahan) / scale, 1e-13)
      << "pairwise=" << pairwise << " kahan=" << kahan;
}

TEST(IntegrateTest, DeterministicAcrossInvocations) {
  // Bit-equality across repeated invocations. At Phase 3 this is trivially
  // satisfied because the reducer is single-threaded; Phase 21 expands this
  // to vary OMP_NUM_THREADS while keeping bit-identical output.
  const int N = 32;
  const double h = 2.0 * kPi / static_cast<double>(N);

  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> f(grid, [&](double x, double y, double z) {
    return std::sin(x) * std::sin(y) * std::sin(z);
  });

  const double a_pair = integrate(f, Pairwise{});
  const double b_pair = integrate(f, Pairwise{});
  const double a_kahan = integrate(f, Kahan{});
  const double b_kahan = integrate(f, Kahan{});

  EXPECT_TRUE(bitEqual(a_pair, b_pair));
  EXPECT_TRUE(bitEqual(a_kahan, b_kahan));
}
