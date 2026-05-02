#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/diff/Laplacian.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
using gridcalc::diff::laplacian;

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

TEST(LaplacianTest, RecoversTrigEigenvalue) {
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  const double kx = 1.0;
  const double ky = 2.0;
  const double kz = 1.0;
  const double eigenvalue = -(kx * kx + ky * ky + kz * kz);  // = -6

  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [&](double x, double y, double z) {
    return std::sin(kx * x) * std::sin(ky * y) * std::sin(kz * z);
  });
  Field<double> expected(grid, [&](double x, double y, double z) {
    return eigenvalue * std::sin(kx * x) * std::sin(ky * y) * std::sin(kz * z);
  });

  Field<double> lap = laplacian(psi);

  // Leading-order relative error for the 2nd-order central-difference
  // eigenvalue on a trig product is (kx^4+ky^4+kz^4) * h^2 / (12 * (kx^2+ky^2+kz^2)).
  // For k=(1,2,1) at h=2pi/32, this evaluates to ~9.6e-3.
  const double rel_err = maxAbsDiff(lap, expected) / maxAbs(expected);
  EXPECT_LT(rel_err, 1.5e-2) << "relative max-norm error was " << rel_err;
}

TEST(LaplacianTest, SecondOrderConvergence) {
  const double L = 2.0 * kPi;
  const double kx = 1.0;
  const double ky = 1.0;
  const double kz = 1.0;
  const double eigenvalue = -(kx * kx + ky * ky + kz * kz);  // = -3

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
    Field<double> lap = laplacian(psi);
    errors.push_back(maxAbsDiff(lap, expected) / maxAbs(expected));
  }

  // Halving h (doubling N) should reduce error by ~4 for a 2nd-order scheme.
  for (std::size_t i = 1; i < errors.size(); ++i) {
    const double ratio = errors[i - 1] / errors[i];
    EXPECT_GT(ratio, 3.5) << "ratio at step " << i << " was " << ratio;
    EXPECT_LT(ratio, 4.5) << "ratio at step " << i << " was " << ratio;
  }

  // Least-squares log-log slope of error vs h. error ~ C h^p, so
  // log(error) = p * log(h) + const, and we expect p in [1.8, 2.2].
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

TEST(LaplacianTest, OutputGridMatchesInput) {
  Grid grid(8, 4, 16, Vec3d(0.1, 0.2, 0.05));
  Field<double> f(grid);
  Field<double> lap = laplacian(f);
  EXPECT_EQ(lap.getGrid().getNx(), grid.getNx());
  EXPECT_EQ(lap.getGrid().getNy(), grid.getNy());
  EXPECT_EQ(lap.getGrid().getNz(), grid.getNz());
  EXPECT_DOUBLE_EQ(lap.getGrid().getCellVolume(), grid.getCellVolume());
}
