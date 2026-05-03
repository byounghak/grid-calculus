#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/diff/Divergence.hpp>
#include <gridcalc/diff/Gradient.hpp>
#include <gridcalc/diff/Laplacian.hpp>
#include <gridcalc/stencil/CentralDifference.hpp>
#include <gridcalc/stencil/FirstDerivative.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
using gridcalc::diff::divergence;
using gridcalc::diff::gradient;
using gridcalc::diff::laplacian;
using gridcalc::stencil::Coefficients;
using gridcalc::stencil::FirstDerivative;

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

double maxComponentAbs(const Field<Vec3d>& f) {
  double m = 0.0;
  for (std::size_t n = 0; n < f.getSize(); ++n) {
    const double a = f.data()[n].cwiseAbs().maxCoeff();
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

}  // namespace

// ---------- weight-table verification ----------

TEST(HigherOrderTest, Coefficients4WeightsAreCorrect) {
  EXPECT_EQ(Coefficients<4>::radius, 2);
  EXPECT_DOUBLE_EQ(Coefficients<4>::weights[0], -1.0 / 12.0);
  EXPECT_DOUBLE_EQ(Coefficients<4>::weights[1], 4.0 / 3.0);
  EXPECT_DOUBLE_EQ(Coefficients<4>::weights[2], -5.0 / 2.0);
  EXPECT_DOUBLE_EQ(Coefficients<4>::weights[3], 4.0 / 3.0);
  EXPECT_DOUBLE_EQ(Coefficients<4>::weights[4], -1.0 / 12.0);

  EXPECT_EQ(FirstDerivative<4>::radius, 2);
  EXPECT_DOUBLE_EQ(FirstDerivative<4>::weights[0], 1.0 / 12.0);
  EXPECT_DOUBLE_EQ(FirstDerivative<4>::weights[1], -2.0 / 3.0);
  EXPECT_DOUBLE_EQ(FirstDerivative<4>::weights[2], 0.0);
  EXPECT_DOUBLE_EQ(FirstDerivative<4>::weights[3], 2.0 / 3.0);
  EXPECT_DOUBLE_EQ(FirstDerivative<4>::weights[4], -1.0 / 12.0);
}

// ---------- order-4 convergence sweeps ----------

TEST(HigherOrderTest, LaplacianOrder4Sweep) {
  const double L = 2.0 * kPi;
  const std::vector<int> Ns = {16, 32, 64};
  std::vector<double> errors;
  for (int N : Ns) {
    const double h = L / static_cast<double>(N);
    Grid grid(N, N, N, Vec3d(h, h, h));
    Field<double> psi(grid, [&](double x, double y, double z) {
      return std::sin(x) * std::sin(y) * std::sin(z);
    });
    Field<double> expected(grid, [&](double x, double y, double z) {
      return -3.0 * std::sin(x) * std::sin(y) * std::sin(z);
    });
    Field<double> lap = laplacian<4>(psi);
    errors.push_back(maxAbsDiff(lap, expected) / maxAbs(expected));
  }
  const double slope = logLogSlope(Ns, errors, L);
  EXPECT_GT(slope, 3.5) << "order-4 Laplacian slope was " << slope;
  EXPECT_LT(slope, 4.5) << "order-4 Laplacian slope was " << slope;
}

TEST(HigherOrderTest, GradientOrder4Sweep) {
  const double L = 2.0 * kPi;
  const std::vector<int> Ns = {16, 32, 64};
  std::vector<double> errors;
  for (int N : Ns) {
    const double h = L / static_cast<double>(N);
    Grid grid(N, N, N, Vec3d(h, h, h));
    Field<double> psi(grid, [&](double x, double y, double z) {
      return std::sin(x) * std::sin(y) * std::sin(z);
    });
    Field<Vec3d> expected(grid, [&](double x, double y, double z) {
      const double sx = std::sin(x), sy = std::sin(y), sz = std::sin(z);
      const double cx = std::cos(x), cy = std::cos(y), cz = std::cos(z);
      return Vec3d(cx * sy * sz, sx * cy * sz, sx * sy * cz);
    });
    Field<Vec3d> grad = gradient<4>(psi);
    errors.push_back(maxComponentAbsDiff(grad, expected) / maxComponentAbs(expected));
  }
  const double slope = logLogSlope(Ns, errors, L);
  EXPECT_GT(slope, 3.5) << "order-4 gradient slope was " << slope;
  EXPECT_LT(slope, 4.5) << "order-4 gradient slope was " << slope;
}

TEST(HigherOrderTest, DivergenceOrder4Sweep) {
  const double L = 2.0 * kPi;
  const std::vector<int> Ns = {16, 32, 64};
  std::vector<double> errors;
  for (int N : Ns) {
    const double h = L / static_cast<double>(N);
    Grid grid(N, N, N, Vec3d(h, h, h));
    Field<Vec3d> V(grid, [&](double x, double y, double z) {
      return Vec3d(std::sin(x) * std::cos(y),
                   std::cos(x) * std::sin(y),
                   std::sin(z));
    });
    Field<double> expected(grid, [&](double x, double y, double z) {
      return 2.0 * std::cos(x) * std::cos(y) + std::cos(z);
    });
    Field<double> div = divergence<4>(V);
    errors.push_back(maxAbsDiff(div, expected) / maxAbs(expected));
  }
  const double slope = logLogSlope(Ns, errors, L);
  EXPECT_GT(slope, 3.5) << "order-4 divergence slope was " << slope;
  EXPECT_LT(slope, 4.5) << "order-4 divergence slope was " << slope;
}

// ---------- order-2 vs order-4 absolute-accuracy comparison ----------

TEST(HigherOrderTest, LaplacianOrder4MoreAccurateThanOrder2) {
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [&](double x, double y, double z) {
    return std::sin(x) * std::sin(y) * std::sin(z);
  });
  Field<double> expected(grid, [&](double x, double y, double z) {
    return -3.0 * std::sin(x) * std::sin(y) * std::sin(z);
  });

  const double scale = maxAbs(expected);
  const double err2 = maxAbsDiff(laplacian<2>(psi), expected) / scale;
  const double err4 = maxAbsDiff(laplacian<4>(psi), expected) / scale;

  EXPECT_LT(err4 * 10.0, err2)
      << "err2=" << err2 << " err4=" << err4 << " (ratio " << err2 / err4 << ")";
}

TEST(HigherOrderTest, GradientOrder4MoreAccurateThanOrder2) {
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [&](double x, double y, double z) {
    return std::sin(x) * std::sin(y) * std::sin(z);
  });
  Field<Vec3d> expected(grid, [&](double x, double y, double z) {
    const double sx = std::sin(x), sy = std::sin(y), sz = std::sin(z);
    const double cx = std::cos(x), cy = std::cos(y), cz = std::cos(z);
    return Vec3d(cx * sy * sz, sx * cy * sz, sx * sy * cz);
  });

  const double scale = maxComponentAbs(expected);
  const double err2 = maxComponentAbsDiff(gradient<2>(psi), expected) / scale;
  const double err4 = maxComponentAbsDiff(gradient<4>(psi), expected) / scale;

  EXPECT_LT(err4 * 10.0, err2)
      << "err2=" << err2 << " err4=" << err4 << " (ratio " << err2 / err4 << ")";
}

TEST(HigherOrderTest, DivergenceOrder4MoreAccurateThanOrder2) {
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<Vec3d> V(grid, [&](double x, double y, double z) {
    return Vec3d(std::sin(x) * std::cos(y),
                 std::cos(x) * std::sin(y),
                 std::sin(z));
  });
  Field<double> expected(grid, [&](double x, double y, double z) {
    return 2.0 * std::cos(x) * std::cos(y) + std::cos(z);
  });

  const double scale = maxAbs(expected);
  const double err2 = maxAbsDiff(divergence<2>(V), expected) / scale;
  const double err4 = maxAbsDiff(divergence<4>(V), expected) / scale;

  EXPECT_LT(err4 * 10.0, err2)
      << "err2=" << err2 << " err4=" << err4 << " (ratio " << err2 / err4 << ")";
}
