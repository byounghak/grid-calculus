// Tests for the spectral wavenumber arrays (kxRfft / kyFull / kzFull).
//
// Pin the numpy.fft.fftfreq convention at every parity. Includes direct
// element-by-element equality against a hand-computed reference for both
// even and odd N, plus a manufactured-solution highest-positive-harmonic
// acceptance test on odd-N grids -- the strongest gate for the 0.14.3
// fix to `kyFull` / `kzFull`.

#ifdef GRIDCALC_USE_FFT

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/spectral/Derivatives.hpp>
#include <gridcalc/spectral/Wavenumbers.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
namespace gs = gridcalc::spectral;

namespace {

constexpr double kPi = 3.14159265358979323846;

// Hand-computed numpy.fftfreq reference, scaled by 2*pi / (N*h):
// indices [0..N-1] map to signed harmonics [0, 1, ..., k_max, -k_max', ..., -1]
// where for even N: k_max = N/2 - 1 and the next slot is the Nyquist -N/2;
// for odd N: k_max = (N-1)/2 and the next slot is -(N-1)/2 (no Nyquist).
std::vector<double> fftfreqReference(int N, double L) {
  const double scale = 2.0 * kPi * static_cast<double>(N) / L /
                       static_cast<double>(N);  // = 2*pi / L
  // (Note: scale = 2*pi / (N*h); with h = L/N, scale = 2*pi/L. The two
  // forms agree.)
  std::vector<double> ref(static_cast<std::size_t>(N));
  for (int n = 0; n < N; ++n) {
    const int n_signed = (2 * n < N) ? n : n - N;
    ref[static_cast<std::size_t>(n)] = scale * static_cast<double>(n_signed);
  }
  return ref;
}

// rfftfreq reference: [0, 1, ..., N/2] all non-negative.
std::vector<double> rfftfreqReference(int N, double L) {
  const double scale = 2.0 * kPi / L;
  std::vector<double> ref(static_cast<std::size_t>(N / 2 + 1));
  for (std::size_t n = 0; n < ref.size(); ++n) {
    ref[n] = scale * static_cast<double>(n);
  }
  return ref;
}

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

}  // namespace

// =============================================================================
// kyFull — direct array equality at even and odd N
// =============================================================================

class KyFullParityTest : public ::testing::TestWithParam<int> {};

TEST_P(KyFullParityTest, MatchesFftfreqReference) {
  const int Ny = GetParam();
  const Grid g = makeUnitDomainGrid(Ny);
  const auto ky = gs::kyFull(g);
  const auto ref = fftfreqReference(Ny, 2.0 * kPi);
  ASSERT_EQ(ky.size(), ref.size());
  for (std::size_t n = 0; n < ky.size(); ++n) {
    EXPECT_DOUBLE_EQ(ky[n], ref[n])
        << "kyFull mismatch at Ny=" << Ny << ", n=" << n;
  }
}

INSTANTIATE_TEST_SUITE_P(EvenN, KyFullParityTest,
                         ::testing::Values(4, 6, 16, 32));
INSTANTIATE_TEST_SUITE_P(OddN, KyFullParityTest,
                         ::testing::Values(5, 7, 15, 31));

// =============================================================================
// kzFull — direct array equality at even and odd N
// =============================================================================

class KzFullParityTest : public ::testing::TestWithParam<int> {};

TEST_P(KzFullParityTest, MatchesFftfreqReference) {
  const int Nz = GetParam();
  const Grid g = makeUnitDomainGrid(Nz);
  const auto kz = gs::kzFull(g);
  const auto ref = fftfreqReference(Nz, 2.0 * kPi);
  ASSERT_EQ(kz.size(), ref.size());
  for (std::size_t n = 0; n < kz.size(); ++n) {
    EXPECT_DOUBLE_EQ(kz[n], ref[n])
        << "kzFull mismatch at Nz=" << Nz << ", n=" << n;
  }
}

INSTANTIATE_TEST_SUITE_P(EvenN, KzFullParityTest,
                         ::testing::Values(4, 6, 16, 32));
INSTANTIATE_TEST_SUITE_P(OddN, KzFullParityTest,
                         ::testing::Values(5, 7, 15, 31));

// =============================================================================
// kxRfft — regression: rfft path is parity-independent
// =============================================================================

class KxRfftParityTest : public ::testing::TestWithParam<int> {};

TEST_P(KxRfftParityTest, MatchesRfftfreqReference) {
  const int Nx = GetParam();
  const Grid g = makeUnitDomainGrid(Nx);
  const auto kx = gs::kxRfft(g);
  const auto ref = rfftfreqReference(Nx, 2.0 * kPi);
  ASSERT_EQ(kx.size(), ref.size());
  for (std::size_t n = 0; n < kx.size(); ++n) {
    EXPECT_DOUBLE_EQ(kx[n], ref[n])
        << "kxRfft mismatch at Nx=" << Nx << ", n=" << n;
  }
}

INSTANTIATE_TEST_SUITE_P(EvenAndOdd, KxRfftParityTest,
                         ::testing::Values(4, 5, 6, 7, 16, 17, 32));

// =============================================================================
// Highest-positive-harmonic acceptance — manufactured-solution gate
// =============================================================================
//
// For an odd-N axis with N samples on a periodic interval [0, L), the
// highest positive harmonic is at index (N-1)/2. Construct a field
//
//   psi(y) = sin(k_max * y),   k_max = 2*pi * ((N-1)/2) / L
//
// applying spectral::laplacian must recover -k_max^2 * psi to round-off
// (sin is an exact eigenfunction of d^2/dy^2 in the discrete Fourier
// basis at this band-limited resolution; spectral differentiation is
// exact by construction).
//
// Without the 0.14.3 fix, kyFull at index (N-1)/2 returns the wrong
// signed wavenumber (e.g., -(N+1)/2 instead of +(N-1)/2 for N=5), so
// |k|^2 at that index is wrong and the recovered Laplacian fails by
// O(|k_max|^2) magnitude. With the fix the test passes to round-off.

namespace {

Field<double> makeYHighestHarmonic(const Grid& g) {
  const int Ny = g.getNy();
  const double L = static_cast<double>(Ny) * g.getCellSize()(1);
  const double k_max = 2.0 * kPi * static_cast<double>((Ny - 1) / 2) / L;
  return Field<double>(g, [k_max](double, double y, double) {
    return std::sin(k_max * y);
  });
}

Field<double> makeZHighestHarmonic(const Grid& g) {
  const int Nz = g.getNz();
  const double L = static_cast<double>(Nz) * g.getCellSize()(2);
  const double k_max = 2.0 * kPi * static_cast<double>((Nz - 1) / 2) / L;
  return Field<double>(g, [k_max](double, double, double z) {
    return std::sin(k_max * z);
  });
}

}  // namespace

class HighestHarmonicAcceptanceY : public ::testing::TestWithParam<int> {};

TEST_P(HighestHarmonicAcceptanceY, SpectralLaplacianRecoversEigenvalue) {
  const int Ny = GetParam();
  const Grid g = makeUnitDomainGrid(Ny);
  const auto psi = makeYHighestHarmonic(g);
  const double L = 2.0 * kPi;
  const double k_max = 2.0 * kPi * static_cast<double>((Ny - 1) / 2) / L;

  // analytical: laplacian(sin(k_max * y)) = -k_max^2 * sin(k_max * y)
  Field<double> expected(g, [k_max](double, double y, double) {
    return -k_max * k_max * std::sin(k_max * y);
  });

  const auto result = gs::laplacian(psi);
  const double err = maxAbsDiff(result, expected);
  const double scale = std::max(maxAbs(expected), 1e-300);
  // Spectral differentiation is exact for band-limited inputs; allow a
  // generous round-off floor scaled by k_max^2 (the magnitude of the
  // eigenvalue applied per cell). At Ny=5, k_max^2 ~ 4; the floor is
  // dominated by the FFT's machine-precision residual.
  EXPECT_LT(err / scale, 1e-12)
      << "Ny=" << Ny << " err=" << err << " k_max^2=" << k_max * k_max;
}

INSTANTIATE_TEST_SUITE_P(OddN, HighestHarmonicAcceptanceY,
                         ::testing::Values(5, 7, 15, 31));

class HighestHarmonicAcceptanceZ : public ::testing::TestWithParam<int> {};

TEST_P(HighestHarmonicAcceptanceZ, SpectralLaplacianRecoversEigenvalue) {
  const int Nz = GetParam();
  const Grid g = makeUnitDomainGrid(Nz);
  const auto psi = makeZHighestHarmonic(g);
  const double L = 2.0 * kPi;
  const double k_max = 2.0 * kPi * static_cast<double>((Nz - 1) / 2) / L;

  Field<double> expected(g, [k_max](double, double, double z) {
    return -k_max * k_max * std::sin(k_max * z);
  });

  const auto result = gs::laplacian(psi);
  const double err = maxAbsDiff(result, expected);
  const double scale = std::max(maxAbs(expected), 1e-300);
  EXPECT_LT(err / scale, 1e-12)
      << "Nz=" << Nz << " err=" << err << " k_max^2=" << k_max * k_max;
}

INSTANTIATE_TEST_SUITE_P(OddN, HighestHarmonicAcceptanceZ,
                         ::testing::Values(5, 7, 15, 31));

#endif  // GRIDCALC_USE_FFT
