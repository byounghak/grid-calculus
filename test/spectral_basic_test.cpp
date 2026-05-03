#ifdef GRIDCALC_USE_FFT

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/spectral/Derivatives.hpp>
#include <gridcalc/spectral/Fft.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;

namespace {

constexpr double kPi = 3.14159265358979323846;

double maxAbsDiff(const Field<double>& a, const Field<double>& b) {
  double m = 0.0;
  for (std::size_t n = 0; n < a.getSize(); ++n) {
    const double d = std::abs(a.data()[n] - b.data()[n]);
    if (d > m) m = d;
  }
  return m;
}

double maxAbs(const Field<double>& f) {
  double m = 0.0;
  for (double v : f) {
    const double a = std::abs(v);
    if (a > m) m = a;
  }
  return m;
}

}  // namespace

TEST(SpectralBasicTest, R2cC2rRoundTripRecoversInput) {
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / N;
  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [](double x, double y, double z) {
    return std::sin(x) * std::sin(y) * std::sin(z);
  });
  auto spectrum = gridcalc::spectral::forwardR2C(psi);
  Field<double> roundtrip = gridcalc::spectral::backwardC2R(spectrum, grid);
  EXPECT_LT(maxAbsDiff(psi, roundtrip), 5e-13);
}

TEST(SpectralBasicTest, LaplacianRecoversNegKSquared) {
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / N;
  Grid grid(N, N, N, Vec3d(h, h, h));
  // psi = sin(x + y + z); k = (1, 1, 1); |k|^2 = 3.
  Field<double> psi(grid, [](double x, double y, double z) {
    return std::sin(x + y + z);
  });
  Field<double> expected(grid, [](double x, double y, double z) {
    return -3.0 * std::sin(x + y + z);
  });
  Field<double> result = gridcalc::spectral::laplacian(psi);
  EXPECT_LT(maxAbsDiff(result, expected) / maxAbs(expected), 5e-12);
}

TEST(SpectralBasicTest, BiharmonicRecoversKToTheFourth) {
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / N;
  Grid grid(N, N, N, Vec3d(h, h, h));
  // psi = sin(x + y + z); |k|^4 = 9.
  Field<double> psi(grid, [](double x, double y, double z) {
    return std::sin(x + y + z);
  });
  Field<double> expected(grid, [](double x, double y, double z) {
    return 9.0 * std::sin(x + y + z);
  });
  Field<double> result = gridcalc::spectral::biharmonic(psi);
  // Tolerance is looser than laplacian (5e-12) because |k|^4 = 9 amplifies
  // round-off vs |k|^2 = 3.
  EXPECT_LT(maxAbsDiff(result, expected) / maxAbs(expected), 5e-11);
}

TEST(SpectralBasicTest, GradientRecoversKCos) {
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / N;
  Grid grid(N, N, N, Vec3d(h, h, h));
  // psi = sin(x + y + z); grad psi = (cos, cos, cos)(x + y + z).
  Field<double> psi(grid, [](double x, double y, double z) {
    return std::sin(x + y + z);
  });
  auto grad = gridcalc::spectral::gradient(psi);
  Field<Vec3d> expected(grid, [](double x, double y, double z) {
    const double c = std::cos(x + y + z);
    return Vec3d(c, c, c);
  });
  double m = 0.0;
  for (std::size_t n = 0; n < grad.getSize(); ++n) {
    const double d = (grad.data()[n] - expected.data()[n]).cwiseAbs().maxCoeff();
    if (d > m) m = d;
  }
  EXPECT_LT(m, 5e-12) << "max component diff = " << m;
}

TEST(SpectralBasicTest, NyquistZeroingForOddRank) {
  // Construct a single-mode field at the x-Nyquist (N = 16, half_N = 8); verify
  // that spectral::partial<3, 0, 0> returns the zero field.
  const int N = 16;
  const double L = 2.0 * kPi;
  const double h = L / N;
  Grid grid(N, N, N, Vec3d(h, h, h));
  // Mode at x-Nyquist: cos((N/2) * x) = cos((Nx/2) * 2*pi/N * i*h) = cos(pi*i)
  // for grid index i — i.e., (-1)^i, a real mode at the Nyquist.
  Field<double> psi(grid, [](double x, double /*y*/, double /*z*/) {
    return std::cos(8.0 * x);  // N/2 with N = 16
  });
  Field<double> result = gridcalc::spectral::partial<3, 0, 0>(psi);
  // After Nyquist zeroing, the spectrum at x-Nyquist is zero, so the
  // entire output is zero.
  EXPECT_LT(maxAbs(result), 5e-13);
}

#endif  // GRIDCALC_USE_FFT
