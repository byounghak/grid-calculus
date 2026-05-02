#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/func/Functional.hpp>
#include <gridcalc/func/Integrate.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
using gridcalc::func::evaluate;
using gridcalc::func::Kahan;
using gridcalc::func::Pairwise;

namespace {

constexpr double kPi = 3.14159265358979323846;

// Hand-computed reference values for psi = sin(x) * sin(y) * sin(z) on [0, 2pi]^3:
//   integral_psi_squared      = pi^3
//   integral_grad_squared     = 3 * pi^3       (sum of three axis-aligned cos^2 sin^2 sin^2 terms)
//   integral_psi_lap_psi      = -3 * pi^3      (since lap psi == -3 * psi)
//   GL with W=1               = (3/2) * pi^3 + (411/64) * pi^3 = (507/64) * pi^3
double pi_cubed() {
  return kPi * kPi * kPi;
}

}  // namespace

TEST(FunctionalTest, GinzburgLandauHandComputed) {
  const int N = 64;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  const double W = 1.0;

  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [&](double x, double y, double z) {
    return std::sin(x) * std::sin(y) * std::sin(z);
  });

  auto gl_density = [W](double p, const Vec3d& g, double) {
    const double dpsi2 = p * p - 1.0;
    return 0.5 * g.squaredNorm() + W * dpsi2 * dpsi2;
  };

  const double F = evaluate(psi, gl_density);
  const double F_ref = (3.0 / 2.0 + 411.0 / 64.0) * pi_cubed();

  const double rel_err = std::abs(F - F_ref) / std::abs(F_ref);
  EXPECT_LT(rel_err, 1e-2) << "F=" << F << " F_ref=" << F_ref << " rel_err=" << rel_err;
}

TEST(FunctionalTest, GinzburgLandauSecondOrderConvergence) {
  const double L = 2.0 * kPi;
  const double W = 1.0;
  const double F_ref = (3.0 / 2.0 + 411.0 / 64.0) * pi_cubed();

  const std::vector<int> Ns = {16, 32, 64};
  std::vector<double> errors;
  errors.reserve(Ns.size());

  for (int N : Ns) {
    const double h = L / static_cast<double>(N);
    Grid grid(N, N, N, Vec3d(h, h, h));
    Field<double> psi(grid, [&](double x, double y, double z) {
      return std::sin(x) * std::sin(y) * std::sin(z);
    });

    auto gl_density = [W](double p, const Vec3d& g, double) {
      const double dpsi2 = p * p - 1.0;
      return 0.5 * g.squaredNorm() + W * dpsi2 * dpsi2;
    };

    const double F = evaluate(psi, gl_density);
    errors.push_back(std::abs(F - F_ref) / std::abs(F_ref));
  }

  for (std::size_t i = 1; i < errors.size(); ++i) {
    const double ratio = errors[i - 1] / errors[i];
    EXPECT_GT(ratio, 3.5) << "ratio at step " << i << " was " << ratio;
    EXPECT_LT(ratio, 4.5) << "ratio at step " << i << " was " << ratio;
  }

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
  const double order = (static_cast<double>(n) * sum_xy - sum_x * sum_y) /
                       (static_cast<double>(n) * sum_xx - sum_x * sum_x);
  EXPECT_GT(order, 1.8) << "log-log convergence order was " << order;
  EXPECT_LT(order, 2.2) << "log-log convergence order was " << order;
}

TEST(FunctionalTest, PsiOnlyArity) {
  // Verifies the 1-arg dispatch path (skips gradient + laplacian).
  // f(psi) = psi^2; integral over [0, 2pi]^3 with psi = sin(x)sin(y)sin(z) is pi^3.
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);

  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [&](double x, double y, double z) {
    return std::sin(x) * std::sin(y) * std::sin(z);
  });

  auto f = [](double p) { return p * p; };
  const double F = evaluate(psi, f);
  EXPECT_NEAR(F, pi_cubed(), 1e-12);
}

TEST(FunctionalTest, PsiAndGradArity) {
  // Verifies the 2-arg dispatch path (skips laplacian).
  // f(psi, grad) = (1/2)|grad|^2; integral is (3/2) pi^3.
  const int N = 64;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);

  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [&](double x, double y, double z) {
    return std::sin(x) * std::sin(y) * std::sin(z);
  });

  auto f = [](double, const Vec3d& g) { return 0.5 * g.squaredNorm(); };
  const double F = evaluate(psi, f);
  const double F_ref = 1.5 * pi_cubed();

  const double rel_err = std::abs(F - F_ref) / std::abs(F_ref);
  EXPECT_LT(rel_err, 1e-2) << "F=" << F << " F_ref=" << F_ref;
}

TEST(FunctionalTest, PsiGradLapArity) {
  // Verifies the 3-arg dispatch path (materializes both gradient and laplacian).
  // f(psi, grad, lap) = psi * lap. For psi = sin(x)sin(y)sin(z), lap psi = -3 psi,
  // so the integrand is -3 * psi^2 with analytical integral -3 * pi^3.
  const int N = 64;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);

  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [&](double x, double y, double z) {
    return std::sin(x) * std::sin(y) * std::sin(z);
  });

  auto f = [](double p, const Vec3d&, double l) { return p * l; };
  const double F = evaluate(psi, f);
  const double F_ref = -3.0 * pi_cubed();

  const double rel_err = std::abs(F - F_ref) / std::abs(F_ref);
  EXPECT_LT(rel_err, 1e-2) << "F=" << F << " F_ref=" << F_ref;
}

TEST(FunctionalTest, KahanTagAgrees) {
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  const double W = 1.0;

  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [&](double x, double y, double z) {
    return std::sin(x) * std::sin(y) * std::sin(z);
  });

  auto gl_density = [W](double p, const Vec3d& g, double) {
    const double dpsi2 = p * p - 1.0;
    return 0.5 * g.squaredNorm() + W * dpsi2 * dpsi2;
  };

  const double F_pair = evaluate(psi, gl_density, Pairwise{});
  const double F_kahan = evaluate(psi, gl_density, Kahan{});
  const double scale = std::max(std::abs(F_pair), std::abs(F_kahan));
  EXPECT_LT(std::abs(F_pair - F_kahan) / scale, 1e-13)
      << "F_pair=" << F_pair << " F_kahan=" << F_kahan;
}
