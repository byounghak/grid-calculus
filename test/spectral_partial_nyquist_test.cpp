// Audit gate: parity-gated Nyquist-zeroing in spectral::partial /
// spectral::gradient on odd-N grids (deferred follow-up from 0.14.3,
// discharged at 0.14.4).
//
// The Nyquist-zeroing path in include/gridcalc/spectral/Derivatives.hpp
// is gated on (N_grid % 2 == 0) AND (per-axis derivative count is odd).
// For odd N there is no Nyquist mode (the highest positive harmonic at
// (N - 1)/2 is a real harmonic, not an aliased Nyquist), so the
// zeroing path must short-circuit and the multiplier must be applied
// fully. The 0.14.3 FD-FFT cross-check at N=17 catches gross
// regressions but does NOT exercise odd-derivative-count operators on
// inputs with content at the highest positive harmonic -- that's the
// gap this file closes.
//
// Each test constructs an input field with content at the highest
// positive harmonic on the axis under test, applies the spectral
// derivative, and compares to the closed-form analytical result.
// Without correct parity-gating these would fail by O(|k_max|^N)
// magnitude; with the current correct code they pass to round-off.

#ifdef GRIDCALC_USE_FFT

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/spectral/Derivatives.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
namespace gs = gridcalc::spectral;

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kRelTol = 1e-12;

Grid makeUnitDomainGrid(int N) {
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  return Grid(N, N, N, Vec3d(h, h, h));
}

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

double kMaxOnAxis(int N) {
  const double L = 2.0 * kPi;
  return 2.0 * kPi * static_cast<double>((N - 1) / 2) / L;
}

}  // namespace

// =============================================================================
// First derivative (odd per-axis count = 1): spectral::partial<1,0,0> etc.
// =============================================================================

class SpectralPartialFirstDerivX : public ::testing::TestWithParam<int> {};

TEST_P(SpectralPartialFirstDerivX, HighestHarmonicRecoversCosine) {
  const int Nx = GetParam();
  const Grid g = makeUnitDomainGrid(Nx);
  const double k = kMaxOnAxis(Nx);
  const Field<double> psi(g, [k](double x, double, double) {
    return std::sin(k * x);
  });
  Field<double> expected(g, [k](double x, double, double) {
    return k * std::cos(k * x);
  });
  const auto result = gs::partial<1, 0, 0>(psi);
  const double err = maxAbsDiff(result, expected);
  const double scale = std::max(maxAbs(expected), 1e-300);
  EXPECT_LT(err / scale, kRelTol)
      << "Nx=" << Nx << " k_max=" << k << " err=" << err;
}

INSTANTIATE_TEST_SUITE_P(OddN, SpectralPartialFirstDerivX,
                         ::testing::Values(5, 7, 15, 31));

class SpectralPartialFirstDerivY : public ::testing::TestWithParam<int> {};

TEST_P(SpectralPartialFirstDerivY, HighestHarmonicRecoversCosine) {
  const int Ny = GetParam();
  const Grid g = makeUnitDomainGrid(Ny);
  const double k = kMaxOnAxis(Ny);
  const Field<double> psi(g, [k](double, double y, double) {
    return std::sin(k * y);
  });
  Field<double> expected(g, [k](double, double y, double) {
    return k * std::cos(k * y);
  });
  const auto result = gs::partial<0, 1, 0>(psi);
  const double err = maxAbsDiff(result, expected);
  const double scale = std::max(maxAbs(expected), 1e-300);
  EXPECT_LT(err / scale, kRelTol)
      << "Ny=" << Ny << " k_max=" << k << " err=" << err;
}

INSTANTIATE_TEST_SUITE_P(OddN, SpectralPartialFirstDerivY,
                         ::testing::Values(5, 7, 15, 31));

class SpectralPartialFirstDerivZ : public ::testing::TestWithParam<int> {};

TEST_P(SpectralPartialFirstDerivZ, HighestHarmonicRecoversCosine) {
  const int Nz = GetParam();
  const Grid g = makeUnitDomainGrid(Nz);
  const double k = kMaxOnAxis(Nz);
  const Field<double> psi(g, [k](double, double, double z) {
    return std::sin(k * z);
  });
  Field<double> expected(g, [k](double, double, double z) {
    return k * std::cos(k * z);
  });
  const auto result = gs::partial<0, 0, 1>(psi);
  const double err = maxAbsDiff(result, expected);
  const double scale = std::max(maxAbs(expected), 1e-300);
  EXPECT_LT(err / scale, kRelTol)
      << "Nz=" << Nz << " k_max=" << k << " err=" << err;
}

INSTANTIATE_TEST_SUITE_P(OddN, SpectralPartialFirstDerivZ,
                         ::testing::Values(5, 7, 15, 31));

// =============================================================================
// Third derivative (odd per-axis count = 3): spectral::partial<3,0,0> etc.
// d^3/dx^3 [sin(k*x)] = -k^3 * cos(k*x)
// =============================================================================

class SpectralPartialThirdDerivX : public ::testing::TestWithParam<int> {};

TEST_P(SpectralPartialThirdDerivX, HighestHarmonicRecoversNegCosine) {
  const int Nx = GetParam();
  const Grid g = makeUnitDomainGrid(Nx);
  const double k = kMaxOnAxis(Nx);
  const Field<double> psi(g, [k](double x, double, double) {
    return std::sin(k * x);
  });
  const double k3 = k * k * k;
  Field<double> expected(g, [k, k3](double x, double, double) {
    return -k3 * std::cos(k * x);
  });
  const auto result = gs::partial<3, 0, 0>(psi);
  const double err = maxAbsDiff(result, expected);
  const double scale = std::max(maxAbs(expected), 1e-300);
  EXPECT_LT(err / scale, kRelTol)
      << "Nx=" << Nx << " k_max^3=" << k3 << " err=" << err;
}

INSTANTIATE_TEST_SUITE_P(OddN, SpectralPartialThirdDerivX,
                         ::testing::Values(5, 7, 15, 31));

class SpectralPartialThirdDerivY : public ::testing::TestWithParam<int> {};

TEST_P(SpectralPartialThirdDerivY, HighestHarmonicRecoversNegCosine) {
  const int Ny = GetParam();
  const Grid g = makeUnitDomainGrid(Ny);
  const double k = kMaxOnAxis(Ny);
  const Field<double> psi(g, [k](double, double y, double) {
    return std::sin(k * y);
  });
  const double k3 = k * k * k;
  Field<double> expected(g, [k, k3](double, double y, double) {
    return -k3 * std::cos(k * y);
  });
  const auto result = gs::partial<0, 3, 0>(psi);
  const double err = maxAbsDiff(result, expected);
  const double scale = std::max(maxAbs(expected), 1e-300);
  EXPECT_LT(err / scale, kRelTol)
      << "Ny=" << Ny << " k_max^3=" << k3 << " err=" << err;
}

INSTANTIATE_TEST_SUITE_P(OddN, SpectralPartialThirdDerivY,
                         ::testing::Values(5, 7, 15, 31));

class SpectralPartialThirdDerivZ : public ::testing::TestWithParam<int> {};

TEST_P(SpectralPartialThirdDerivZ, HighestHarmonicRecoversNegCosine) {
  const int Nz = GetParam();
  const Grid g = makeUnitDomainGrid(Nz);
  const double k = kMaxOnAxis(Nz);
  const Field<double> psi(g, [k](double, double, double z) {
    return std::sin(k * z);
  });
  const double k3 = k * k * k;
  Field<double> expected(g, [k, k3](double, double, double z) {
    return -k3 * std::cos(k * z);
  });
  const auto result = gs::partial<0, 0, 3>(psi);
  const double err = maxAbsDiff(result, expected);
  const double scale = std::max(maxAbs(expected), 1e-300);
  EXPECT_LT(err / scale, kRelTol)
      << "Nz=" << Nz << " k_max^3=" << k3 << " err=" << err;
}

INSTANTIATE_TEST_SUITE_P(OddN, SpectralPartialThirdDerivZ,
                         ::testing::Values(5, 7, 15, 31));

// =============================================================================
// spectral::gradient on odd-N grids with content at the highest harmonic
// =============================================================================
//
// gradient returns a Field<Vec3d>; the active-axis component must match
// the closed-form first derivative. The other two components must be
// zero to round-off (no spatial variation along their axes).

namespace {

double maxComponent(const Field<Vec3d>& f, int comp) {
  double m = 0.0;
  for (std::size_t n = 0; n < f.getSize(); ++n) {
    const double v = std::abs(f.data()[n](comp));
    if (v > m) m = v;
  }
  return m;
}

double maxComponentDiff(const Field<Vec3d>& f, int comp,
                        const Field<double>& expected) {
  double m = 0.0;
  for (std::size_t n = 0; n < f.getSize(); ++n) {
    const double d = std::abs(f.data()[n](comp) - expected.data()[n]);
    if (d > m) m = d;
  }
  return m;
}

}  // namespace

class SpectralGradientHighestHarmonicX : public ::testing::TestWithParam<int> {};

TEST_P(SpectralGradientHighestHarmonicX, ActiveComponentRecoversAnalytical) {
  const int Nx = GetParam();
  const Grid g = makeUnitDomainGrid(Nx);
  const double k = kMaxOnAxis(Nx);
  const Field<double> psi(g, [k](double x, double, double) {
    return std::sin(k * x);
  });
  Field<double> expected(g, [k](double x, double, double) {
    return k * std::cos(k * x);
  });
  const auto grad = gs::gradient(psi);
  const double err_x = maxComponentDiff(grad, 0, expected);
  const double scale_x = std::max(maxAbs(expected), 1e-300);
  EXPECT_LT(err_x / scale_x, kRelTol)
      << "x-component err at Nx=" << Nx << ": " << err_x;
  // y and z components must be zero (no spatial variation on those axes)
  EXPECT_LT(maxComponent(grad, 1), 1e-10)
      << "y-component should be ~0 at Nx=" << Nx;
  EXPECT_LT(maxComponent(grad, 2), 1e-10)
      << "z-component should be ~0 at Nx=" << Nx;
}

INSTANTIATE_TEST_SUITE_P(OddN, SpectralGradientHighestHarmonicX,
                         ::testing::Values(5, 7, 15, 31));

class SpectralGradientHighestHarmonicY : public ::testing::TestWithParam<int> {};

TEST_P(SpectralGradientHighestHarmonicY, ActiveComponentRecoversAnalytical) {
  const int Ny = GetParam();
  const Grid g = makeUnitDomainGrid(Ny);
  const double k = kMaxOnAxis(Ny);
  const Field<double> psi(g, [k](double, double y, double) {
    return std::sin(k * y);
  });
  Field<double> expected(g, [k](double, double y, double) {
    return k * std::cos(k * y);
  });
  const auto grad = gs::gradient(psi);
  const double err_y = maxComponentDiff(grad, 1, expected);
  const double scale_y = std::max(maxAbs(expected), 1e-300);
  EXPECT_LT(err_y / scale_y, kRelTol)
      << "y-component err at Ny=" << Ny << ": " << err_y;
  EXPECT_LT(maxComponent(grad, 0), 1e-10)
      << "x-component should be ~0 at Ny=" << Ny;
  EXPECT_LT(maxComponent(grad, 2), 1e-10)
      << "z-component should be ~0 at Ny=" << Ny;
}

INSTANTIATE_TEST_SUITE_P(OddN, SpectralGradientHighestHarmonicY,
                         ::testing::Values(5, 7, 15, 31));

class SpectralGradientHighestHarmonicZ : public ::testing::TestWithParam<int> {};

TEST_P(SpectralGradientHighestHarmonicZ, ActiveComponentRecoversAnalytical) {
  const int Nz = GetParam();
  const Grid g = makeUnitDomainGrid(Nz);
  const double k = kMaxOnAxis(Nz);
  const Field<double> psi(g, [k](double, double, double z) {
    return std::sin(k * z);
  });
  Field<double> expected(g, [k](double, double, double z) {
    return k * std::cos(k * z);
  });
  const auto grad = gs::gradient(psi);
  const double err_z = maxComponentDiff(grad, 2, expected);
  const double scale_z = std::max(maxAbs(expected), 1e-300);
  EXPECT_LT(err_z / scale_z, kRelTol)
      << "z-component err at Nz=" << Nz << ": " << err_z;
  EXPECT_LT(maxComponent(grad, 0), 1e-10)
      << "x-component should be ~0 at Nz=" << Nz;
  EXPECT_LT(maxComponent(grad, 1), 1e-10)
      << "y-component should be ~0 at Nz=" << Nz;
}

INSTANTIATE_TEST_SUITE_P(OddN, SpectralGradientHighestHarmonicZ,
                         ::testing::Values(5, 7, 15, 31));

#endif  // GRIDCALC_USE_FFT
