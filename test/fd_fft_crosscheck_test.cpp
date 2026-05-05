// FD-FFT cross-check fixture (Phase 9). For every Phase 1-8 finite-difference
// operator at every advertised accuracy order, this fixture verifies that the
// FD-FFT discrepancy on a smooth band-limited input scales as O(h^p) where p
// is the FD order. Becomes a permanent CI gate: future FD changes that
// silently break accuracy fail this test.

#ifdef GRIDCALC_USE_FFT

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/diff/Biharmonic.hpp>
#include <gridcalc/diff/Divergence.hpp>
#include <gridcalc/diff/FourthOrder.hpp>
#include <gridcalc/diff/Gradient.hpp>
#include <gridcalc/diff/Laplacian.hpp>
#include <gridcalc/diff/MixedPartial.hpp>
#include <gridcalc/diff/ThirdOrder.hpp>
#include <gridcalc/spectral/Derivatives.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
namespace gd = gridcalc::diff;
namespace gs = gridcalc::spectral;

namespace {

constexpr double kPi = 3.14159265358979323846;
// Sweep includes one odd-N grid (N=17) starting at 0.14.3, so every
// FD<->FFT slope-band assertion is exercised on at least one odd-N
// resolution. h=2*pi/17 ~ 0.37 keeps the truncation-analysis regime
// intact (unlike N=5 at h~1.257, which would invalidate the slope
// test -- see specs/2026-05-04-fix-stencil-aliasing-on-small-axes/).
const std::vector<int> kSweepNs = {16, 17, 32, 64, 128};

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

double logLogSlope(const std::vector<int>& Ns, const std::vector<double>& errs,
                   double L) {
  const std::size_t n = Ns.size();
  double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double h = L / static_cast<double>(Ns[i]);
    const double x = std::log(h);
    const double y = std::log(errs[i]);
    sum_x += x;
    sum_y += y;
    sum_xy += x * y;
    sum_xx += x * x;
  }
  return (static_cast<double>(n) * sum_xy - sum_x * sum_y) /
         (static_cast<double>(n) * sum_xx - sum_x * sum_x);
}

Field<double> manufacturedPsi(const Grid& grid) {
  return Field<double>(grid, [](double x, double y, double z) {
    return std::sin(x + y + z);
  });
}

template <int Nx, int Ny, int Nz, typename FdOp>
double crossCheckSlopePartial(FdOp&& fd_op) {
  const double L = 2.0 * kPi;
  std::vector<double> errs;
  errs.reserve(kSweepNs.size());
  for (int N : kSweepNs) {
    const double h = L / static_cast<double>(N);
    Grid grid(N, N, N, Vec3d(h, h, h));
    Field<double> psi = manufacturedPsi(grid);
    Field<double> fd_result = fd_op(psi);
    Field<double> sp_result = gs::partial<Nx, Ny, Nz>(psi);
    errs.push_back(maxAbsDiff(fd_result, sp_result) /
                   std::max(maxAbs(sp_result), 1e-300));
  }
  return logLogSlope(kSweepNs, errs, L);
}

template <typename FdOp, typename SpOp>
double crossCheckSlopeNamed(FdOp&& fd_op, SpOp&& sp_op) {
  const double L = 2.0 * kPi;
  std::vector<double> errs;
  errs.reserve(kSweepNs.size());
  for (int N : kSweepNs) {
    const double h = L / static_cast<double>(N);
    Grid grid(N, N, N, Vec3d(h, h, h));
    Field<double> psi = manufacturedPsi(grid);
    Field<double> fd_result = fd_op(psi);
    Field<double> sp_result = sp_op(psi);
    errs.push_back(maxAbsDiff(fd_result, sp_result) /
                   std::max(maxAbs(sp_result), 1e-300));
  }
  return logLogSlope(kSweepNs, errs, L);
}

double maxComponentAbsDiff(const Field<Vec3d>& a, const Field<Vec3d>& b) {
  double m = 0.0;
  for (std::size_t n = 0; n < a.getSize(); ++n) {
    const double d = (a.data()[n] - b.data()[n]).cwiseAbs().maxCoeff();
    if (d > m) m = d;
  }
  return m;
}

double maxComponentAbs(const Field<Vec3d>& f) {
  double m = 0.0;
  for (std::size_t n = 0; n < f.getSize(); ++n) {
    const double a = f.data()[n].cwiseAbs().maxCoeff();
    if (a > m) m = a;
  }
  return m;
}

template <typename GradOp>
double crossCheckSlopeGradient(GradOp&& grad_op) {
  const double L = 2.0 * kPi;
  std::vector<double> errs;
  errs.reserve(kSweepNs.size());
  for (int N : kSweepNs) {
    const double h = L / static_cast<double>(N);
    Grid grid(N, N, N, Vec3d(h, h, h));
    Field<double> psi = manufacturedPsi(grid);
    Field<Vec3d> fd_result = grad_op(psi);
    Field<Vec3d> sp_result = gs::gradient(psi);
    errs.push_back(maxComponentAbsDiff(fd_result, sp_result) /
                   std::max(maxComponentAbs(sp_result), 1e-300));
  }
  return logLogSlope(kSweepNs, errs, L);
}

template <typename DivOp>
double crossCheckSlopeDivergence(DivOp&& div_op) {
  const double L = 2.0 * kPi;
  std::vector<double> errs;
  errs.reserve(kSweepNs.size());
  for (int N : kSweepNs) {
    const double h = L / static_cast<double>(N);
    Grid grid(N, N, N, Vec3d(h, h, h));
    // Vector field V = (sin(x+y+z), sin(x+y+z), sin(x+y+z));
    // div V = 3 cos(x+y+z). Spectral comparison via the sum of three
    // spectral first-derivatives applied per component.
    Field<Vec3d> V(grid, [](double x, double y, double z) {
      const double s = std::sin(x + y + z);
      return Vec3d(s, s, s);
    });
    Field<double> fd_result = div_op(V);

    // Build the spectral reference: scalar first-derivatives of each
    // component, summed. V's components are identical here (V = (s, s, s)),
    // so spectral div = partial<1,0,0>(s) + partial<0,1,0>(s) + partial<0,0,1>(s).
    Field<double> psi(grid, [](double x, double y, double z) {
      return std::sin(x + y + z);
    });
    Field<double> dx = gs::partial<1, 0, 0>(psi);
    Field<double> dy = gs::partial<0, 1, 0>(psi);
    Field<double> dz = gs::partial<0, 0, 1>(psi);
    Field<double> sp_result(grid);
    for (std::size_t n = 0; n < sp_result.getSize(); ++n) {
      sp_result.data()[n] = dx.data()[n] + dy.data()[n] + dz.data()[n];
    }
    errs.push_back(maxAbsDiff(fd_result, sp_result) /
                   std::max(maxAbs(sp_result), 1e-300));
  }
  return logLogSlope(kSweepNs, errs, L);
}

}  // namespace

// =============================================================================
// Phase 1: laplacian — vs spectral::laplacian.
// =============================================================================

#define EXPECT_FD_FFT_NAMED(TAG, FD_EXPR, SP_EXPR, ORDER)                       \
  TEST(FdFftCrossCheck, TAG##_Order##ORDER) {                                   \
    const double slope = crossCheckSlopeNamed(                                   \
        [](const Field<double>& f) { return FD_EXPR; },                          \
        [](const Field<double>& f) { return SP_EXPR; });                         \
    EXPECT_GT(slope, ORDER - 0.5) << #TAG "<" #ORDER "> slope=" << slope;        \
    EXPECT_LT(slope, ORDER + 0.5) << #TAG "<" #ORDER "> slope=" << slope;        \
  }

EXPECT_FD_FFT_NAMED(Laplacian, gd::laplacian<2>(f), gs::laplacian(f), 2)
EXPECT_FD_FFT_NAMED(Laplacian, gd::laplacian<4>(f), gs::laplacian(f), 4)

// =============================================================================
// Phase 2: gradient (Field<Vec3d>) and divergence.
// =============================================================================

#define EXPECT_FD_FFT_GRADIENT(ORDER)                                            \
  TEST(FdFftCrossCheck, Gradient_Order##ORDER) {                                 \
    const double slope = crossCheckSlopeGradient(                                \
        [](const Field<double>& f) { return gd::gradient<ORDER>(f); });          \
    EXPECT_GT(slope, ORDER - 0.5) << "gradient<" #ORDER "> slope=" << slope;     \
    EXPECT_LT(slope, ORDER + 0.5) << "gradient<" #ORDER "> slope=" << slope;     \
  }

#define EXPECT_FD_FFT_DIVERGENCE(ORDER)                                          \
  TEST(FdFftCrossCheck, Divergence_Order##ORDER) {                               \
    const double slope = crossCheckSlopeDivergence(                              \
        [](const Field<Vec3d>& v) { return gd::divergence<ORDER>(v); });         \
    EXPECT_GT(slope, ORDER - 0.5) << "divergence<" #ORDER "> slope=" << slope;   \
    EXPECT_LT(slope, ORDER + 0.5) << "divergence<" #ORDER "> slope=" << slope;   \
  }

EXPECT_FD_FFT_GRADIENT(2)
EXPECT_FD_FFT_GRADIENT(4)
EXPECT_FD_FFT_DIVERGENCE(2)
EXPECT_FD_FFT_DIVERGENCE(4)

// =============================================================================
// Phase 8: 31 d-prefix partials × 2 orders + biharmonic × 2 orders.
// =============================================================================

#define EXPECT_FD_FFT_PARTIAL(NAME, NX, NY, NZ, ORDER)                           \
  TEST(FdFftCrossCheck, NAME##_Order##ORDER) {                                   \
    const double slope = crossCheckSlopePartial<NX, NY, NZ>(                     \
        [](const Field<double>& f) { return gd::NAME<ORDER>(f); });              \
    EXPECT_GT(slope, ORDER - 0.5) << #NAME "<" #ORDER "> slope=" << slope;       \
    EXPECT_LT(slope, ORDER + 0.5) << #NAME "<" #ORDER "> slope=" << slope;       \
  }

// Rank-2 partials.
EXPECT_FD_FFT_PARTIAL(dxx, 2, 0, 0, 2)
EXPECT_FD_FFT_PARTIAL(dxx, 2, 0, 0, 4)
EXPECT_FD_FFT_PARTIAL(dyy, 0, 2, 0, 2)
EXPECT_FD_FFT_PARTIAL(dyy, 0, 2, 0, 4)
EXPECT_FD_FFT_PARTIAL(dzz, 0, 0, 2, 2)
EXPECT_FD_FFT_PARTIAL(dzz, 0, 0, 2, 4)
EXPECT_FD_FFT_PARTIAL(dxy, 1, 1, 0, 2)
EXPECT_FD_FFT_PARTIAL(dxy, 1, 1, 0, 4)
EXPECT_FD_FFT_PARTIAL(dxz, 1, 0, 1, 2)
EXPECT_FD_FFT_PARTIAL(dxz, 1, 0, 1, 4)
EXPECT_FD_FFT_PARTIAL(dyz, 0, 1, 1, 2)
EXPECT_FD_FFT_PARTIAL(dyz, 0, 1, 1, 4)

// Rank-3 partials.
EXPECT_FD_FFT_PARTIAL(d3dx3, 3, 0, 0, 2)
EXPECT_FD_FFT_PARTIAL(d3dx3, 3, 0, 0, 4)
EXPECT_FD_FFT_PARTIAL(d3dy3, 0, 3, 0, 2)
EXPECT_FD_FFT_PARTIAL(d3dy3, 0, 3, 0, 4)
EXPECT_FD_FFT_PARTIAL(d3dz3, 0, 0, 3, 2)
EXPECT_FD_FFT_PARTIAL(d3dz3, 0, 0, 3, 4)
EXPECT_FD_FFT_PARTIAL(d3dx2dy, 2, 1, 0, 2)
EXPECT_FD_FFT_PARTIAL(d3dx2dy, 2, 1, 0, 4)
EXPECT_FD_FFT_PARTIAL(d3dxdy2, 1, 2, 0, 2)
EXPECT_FD_FFT_PARTIAL(d3dxdy2, 1, 2, 0, 4)
EXPECT_FD_FFT_PARTIAL(d3dx2dz, 2, 0, 1, 2)
EXPECT_FD_FFT_PARTIAL(d3dx2dz, 2, 0, 1, 4)
EXPECT_FD_FFT_PARTIAL(d3dxdz2, 1, 0, 2, 2)
EXPECT_FD_FFT_PARTIAL(d3dxdz2, 1, 0, 2, 4)
EXPECT_FD_FFT_PARTIAL(d3dy2dz, 0, 2, 1, 2)
EXPECT_FD_FFT_PARTIAL(d3dy2dz, 0, 2, 1, 4)
EXPECT_FD_FFT_PARTIAL(d3dydz2, 0, 1, 2, 2)
EXPECT_FD_FFT_PARTIAL(d3dydz2, 0, 1, 2, 4)
EXPECT_FD_FFT_PARTIAL(d3dxdydz, 1, 1, 1, 2)
EXPECT_FD_FFT_PARTIAL(d3dxdydz, 1, 1, 1, 4)

// Rank-4 partials.
EXPECT_FD_FFT_PARTIAL(d4dx4, 4, 0, 0, 2)
EXPECT_FD_FFT_PARTIAL(d4dx4, 4, 0, 0, 4)
EXPECT_FD_FFT_PARTIAL(d4dy4, 0, 4, 0, 2)
EXPECT_FD_FFT_PARTIAL(d4dy4, 0, 4, 0, 4)
EXPECT_FD_FFT_PARTIAL(d4dz4, 0, 0, 4, 2)
EXPECT_FD_FFT_PARTIAL(d4dz4, 0, 0, 4, 4)
EXPECT_FD_FFT_PARTIAL(d4dx3dy, 3, 1, 0, 2)
EXPECT_FD_FFT_PARTIAL(d4dx3dy, 3, 1, 0, 4)
EXPECT_FD_FFT_PARTIAL(d4dxdy3, 1, 3, 0, 2)
EXPECT_FD_FFT_PARTIAL(d4dxdy3, 1, 3, 0, 4)
EXPECT_FD_FFT_PARTIAL(d4dx3dz, 3, 0, 1, 2)
EXPECT_FD_FFT_PARTIAL(d4dx3dz, 3, 0, 1, 4)
EXPECT_FD_FFT_PARTIAL(d4dxdz3, 1, 0, 3, 2)
EXPECT_FD_FFT_PARTIAL(d4dxdz3, 1, 0, 3, 4)
EXPECT_FD_FFT_PARTIAL(d4dy3dz, 0, 3, 1, 2)
EXPECT_FD_FFT_PARTIAL(d4dy3dz, 0, 3, 1, 4)
EXPECT_FD_FFT_PARTIAL(d4dydz3, 0, 1, 3, 2)
EXPECT_FD_FFT_PARTIAL(d4dydz3, 0, 1, 3, 4)
EXPECT_FD_FFT_PARTIAL(d4dx2dy2, 2, 2, 0, 2)
EXPECT_FD_FFT_PARTIAL(d4dx2dy2, 2, 2, 0, 4)
EXPECT_FD_FFT_PARTIAL(d4dx2dz2, 2, 0, 2, 2)
EXPECT_FD_FFT_PARTIAL(d4dx2dz2, 2, 0, 2, 4)
EXPECT_FD_FFT_PARTIAL(d4dy2dz2, 0, 2, 2, 2)
EXPECT_FD_FFT_PARTIAL(d4dy2dz2, 0, 2, 2, 4)
EXPECT_FD_FFT_PARTIAL(d4dx2dydz, 2, 1, 1, 2)
EXPECT_FD_FFT_PARTIAL(d4dx2dydz, 2, 1, 1, 4)
EXPECT_FD_FFT_PARTIAL(d4dxdy2dz, 1, 2, 1, 2)
EXPECT_FD_FFT_PARTIAL(d4dxdy2dz, 1, 2, 1, 4)
EXPECT_FD_FFT_PARTIAL(d4dxdydz2, 1, 1, 2, 2)
EXPECT_FD_FFT_PARTIAL(d4dxdydz2, 1, 1, 2, 4)

// Phase 8: biharmonic.
EXPECT_FD_FFT_NAMED(Biharmonic, gd::biharmonic<2>(f), gs::biharmonic(f), 2)
EXPECT_FD_FFT_NAMED(Biharmonic, gd::biharmonic<4>(f), gs::biharmonic(f), 4)

#endif  // GRIDCALC_USE_FFT
