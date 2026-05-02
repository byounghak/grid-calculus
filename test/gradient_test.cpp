#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/diff/Gradient.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
using gridcalc::diff::gradient;

namespace {

constexpr double kPi = 3.14159265358979323846;

double maxComponentAbs(const Field<Vec3d>& f) {
  double m = 0.0;
  for (std::size_t n = 0; n < f.getSize(); ++n) {
    const Vec3d& v = f.data()[n];
    const double a = v.cwiseAbs().maxCoeff();
    if (a > m) m = a;
  }
  return m;
}

double maxComponentAbsDiff(const Field<Vec3d>& a, const Field<Vec3d>& b) {
  double m = 0.0;
  for (std::size_t n = 0; n < a.getSize(); ++n) {
    const double d = (a.data()[n] - b.data()[n]).cwiseAbs().maxCoeff();
    if (d > m) m = d;
  }
  return m;
}

}  // namespace

TEST(GradientTest, RecoversTrigGradient) {
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  const double kx = 1.0;
  const double ky = 1.0;
  const double kz = 1.0;

  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [&](double x, double y, double z) {
    return std::sin(kx * x) * std::sin(ky * y) * std::sin(kz * z);
  });
  Field<Vec3d> expected(grid, [&](double x, double y, double z) {
    const double sx = std::sin(kx * x);
    const double sy = std::sin(ky * y);
    const double sz = std::sin(kz * z);
    const double cx = std::cos(kx * x);
    const double cy = std::cos(ky * y);
    const double cz = std::cos(kz * z);
    return Vec3d(kx * cx * sy * sz, ky * sx * cy * sz, kz * sx * sy * cz);
  });

  Field<Vec3d> grad = gradient(psi);

  // Leading-order relative max-norm error for the 2nd-order central-difference
  // gradient on a unit-wavenumber trig product is k^2 * h^2 / 6.
  // For k=1 at h=2*pi/32, this evaluates to ~6.4e-3.
  const double rel_err = maxComponentAbsDiff(grad, expected) / maxComponentAbs(expected);
  EXPECT_LT(rel_err, 1.5e-2) << "relative max-norm error was " << rel_err;
}

TEST(GradientTest, SecondOrderConvergence) {
  const double L = 2.0 * kPi;
  const double kx = 1.0;
  const double ky = 1.0;
  const double kz = 1.0;

  const std::vector<int> Ns = {16, 32, 64};
  std::vector<double> errors;
  errors.reserve(Ns.size());

  for (int N : Ns) {
    const double h = L / static_cast<double>(N);
    Grid grid(N, N, N, Vec3d(h, h, h));
    Field<double> psi(grid, [&](double x, double y, double z) {
      return std::sin(kx * x) * std::sin(ky * y) * std::sin(kz * z);
    });
    Field<Vec3d> expected(grid, [&](double x, double y, double z) {
      const double sx = std::sin(kx * x);
      const double sy = std::sin(ky * y);
      const double sz = std::sin(kz * z);
      const double cx = std::cos(kx * x);
      const double cy = std::cos(ky * y);
      const double cz = std::cos(kz * z);
      return Vec3d(kx * cx * sy * sz, ky * sx * cy * sz, kz * sx * sy * cz);
    });
    Field<Vec3d> grad = gradient(psi);
    errors.push_back(maxComponentAbsDiff(grad, expected) / maxComponentAbs(expected));
  }

  // Halving h (doubling N) should reduce error by ~4 for a 2nd-order scheme.
  for (std::size_t i = 1; i < errors.size(); ++i) {
    const double ratio = errors[i - 1] / errors[i];
    EXPECT_GT(ratio, 3.5) << "ratio at step " << i << " was " << ratio;
    EXPECT_LT(ratio, 4.5) << "ratio at step " << i << " was " << ratio;
  }

  // Least-squares log-log slope of error vs h.
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

TEST(GradientTest, OutputGridMatchesInput) {
  Grid grid(8, 4, 16, Vec3d(0.1, 0.2, 0.05));
  Field<double> f(grid);
  Field<Vec3d> grad = gradient(f);
  EXPECT_EQ(grad.getGrid().getNx(), grid.getNx());
  EXPECT_EQ(grad.getGrid().getNy(), grid.getNy());
  EXPECT_EQ(grad.getGrid().getNz(), grid.getNz());
  EXPECT_DOUBLE_EQ(grad.getGrid().getCellVolume(), grid.getCellVolume());
}
