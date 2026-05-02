#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/diff/Divergence.hpp>
#include <gridcalc/diff/Gradient.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
using gridcalc::diff::divergence;
using gridcalc::diff::gradient;

namespace {

constexpr double kPi = 3.14159265358979323846;

double maxAbs(const Field<double>& f) {
  double m = 0.0;
  for (double v : f) {
    const double av = std::abs(v);
    if (av > m) m = av;
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

}  // namespace

TEST(DivergenceTest, RecoversAnalyticalDivergence) {
  // V = (sin x cos y, cos x sin y, sin z) on the periodic 2pi cube.
  // div V = 2 cos(x) cos(y) + cos(z).
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);

  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<Vec3d> V(grid, [&](double x, double y, double z) {
    return Vec3d(std::sin(x) * std::cos(y), std::cos(x) * std::sin(y), std::sin(z));
  });
  Field<double> expected(grid, [&](double x, double y, double z) {
    return 2.0 * std::cos(x) * std::cos(y) + std::cos(z);
  });

  Field<double> div = divergence(V);

  // Leading-order relative max-norm error is ~h^2/6 per axis at unit
  // wavenumber. For h = 2*pi/32 this is ~6.4e-3; tolerate up to 1.5e-2.
  const double rel_err = maxAbsDiff(div, expected) / maxAbs(expected);
  EXPECT_LT(rel_err, 1.5e-2) << "relative max-norm error was " << rel_err;
}

TEST(DivergenceTest, SecondOrderConvergence) {
  const double L = 2.0 * kPi;
  const std::vector<int> Ns = {16, 32, 64};
  std::vector<double> errors;
  errors.reserve(Ns.size());

  for (int N : Ns) {
    const double h = L / static_cast<double>(N);
    Grid grid(N, N, N, Vec3d(h, h, h));
    Field<Vec3d> V(grid, [&](double x, double y, double z) {
      return Vec3d(std::sin(x) * std::cos(y), std::cos(x) * std::sin(y), std::sin(z));
    });
    Field<double> expected(grid, [&](double x, double y, double z) {
      return 2.0 * std::cos(x) * std::cos(y) + std::cos(z);
    });
    Field<double> div = divergence(V);
    errors.push_back(maxAbsDiff(div, expected) / maxAbs(expected));
  }

  for (std::size_t i = 1; i < errors.size(); ++i) {
    const double ratio = errors[i - 1] / errors[i];
    EXPECT_GT(ratio, 3.5) << "ratio at step " << i << " was " << ratio;
    EXPECT_LT(ratio, 4.5) << "ratio at step " << i << " was " << ratio;
  }

  const std::size_t n = Ns.size();
  double sum_x = 0.0;
  double sum_y = 0.0;
  double sum_xy = 0.0;
  double sum_xx = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double h = L / static_cast<double>(Ns[i]);
    const double x = std::log(h);
    const double y = std::log(errors[i]);
    sum_x += x;
    sum_y += y;
    sum_xy += x * y;
    sum_xx += x * x;
  }
  const double order = (static_cast<double>(n) * sum_xy - sum_x * sum_y) /
                       (static_cast<double>(n) * sum_xx - sum_x * sum_x);
  EXPECT_GT(order, 1.8) << "log-log convergence order was " << order;
  EXPECT_LT(order, 2.2) << "log-log convergence order was " << order;
}

TEST(DivergenceTest, RoundTripWithGradientMatchesLaplacian) {
  // div(grad(psi)) and the analytical Laplacian agree at order 2 over h.
  // Note that D1 o D1 != D2 at finite h (D1 o D1 uses a stride-2 stencil),
  // so we compare to the analytical Laplacian, not to laplacian(psi).
  const double L = 2.0 * kPi;
  const double kx = 1.0;
  const double ky = 1.0;
  const double kz = 1.0;
  const double eigenvalue = -(kx * kx + ky * ky + kz * kz);

  const std::vector<int> Ns = {16, 32, 64};
  std::vector<double> errors;
  errors.reserve(Ns.size());

  for (int N : Ns) {
    const double h = L / static_cast<double>(N);
    Grid grid(N, N, N, Vec3d(h, h, h));
    Field<double> psi(grid, [&](double x, double y, double z) {
      return std::sin(kx * x) * std::sin(ky * y) * std::sin(kz * z);
    });
    Field<double> expected(grid, [&](double x, double y, double z) {
      return eigenvalue * std::sin(kx * x) * std::sin(ky * y) * std::sin(kz * z);
    });
    Field<double> lap_indirect = divergence(gradient(psi));
    errors.push_back(maxAbsDiff(lap_indirect, expected) / maxAbs(expected));
  }

  for (std::size_t i = 1; i < errors.size(); ++i) {
    const double ratio = errors[i - 1] / errors[i];
    EXPECT_GT(ratio, 3.5) << "ratio at step " << i << " was " << ratio;
    EXPECT_LT(ratio, 4.5) << "ratio at step " << i << " was " << ratio;
  }

  const std::size_t n = Ns.size();
  double sum_x = 0.0;
  double sum_y = 0.0;
  double sum_xy = 0.0;
  double sum_xx = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double h = L / static_cast<double>(Ns[i]);
    const double x = std::log(h);
    const double y = std::log(errors[i]);
    sum_x += x;
    sum_y += y;
    sum_xy += x * y;
    sum_xx += x * x;
  }
  const double order = (static_cast<double>(n) * sum_xy - sum_x * sum_y) /
                       (static_cast<double>(n) * sum_xx - sum_x * sum_x);
  EXPECT_GT(order, 1.8) << "log-log convergence order was " << order;
  EXPECT_LT(order, 2.2) << "log-log convergence order was " << order;
}

TEST(DivergenceTest, OutputGridMatchesInput) {
  Grid grid(8, 4, 16, Vec3d(0.1, 0.2, 0.05));
  // Broadcast-init to a zero Vec3d so divergence reads defined values; Eigen's
  // default ctor leaves coefficients uninitialized in release builds.
  Field<Vec3d> V(grid, Vec3d::Zero());
  Field<double> div = divergence(V);
  EXPECT_EQ(div.getGrid().getNx(), grid.getNx());
  EXPECT_EQ(div.getGrid().getNy(), grid.getNy());
  EXPECT_EQ(div.getGrid().getNz(), grid.getNz());
  EXPECT_DOUBLE_EQ(div.getGrid().getCellVolume(), grid.getCellVolume());
}
