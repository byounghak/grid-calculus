#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <vector>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/solve/Diffusion.hpp>
#include <gridcalc/solve/Integrator.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
using gridcalc::solve::diffuse;
using gridcalc::solve::ExplicitEuler;
using gridcalc::solve::explicitEuler;
using gridcalc::solve::integrate;
using gridcalc::solve::RK4;

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

bool bitEqualField(const Field<double>& a, const Field<double>& b) {
  if (a.getSize() != b.getSize()) return false;
  return std::memcmp(a.data(), b.data(), a.getSize() * sizeof(double)) == 0;
}

}  // namespace

TEST(IntegratorTest, RK4MatchesAnalyticAtFixedN) {
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  const double D = 0.1;
  const double dt = 0.05;
  const int n_steps = 100;
  const double T = static_cast<double>(n_steps) * dt;  // 5.0

  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [&](double x, double y, double z) {
    return std::sin(x) * std::sin(y) * std::sin(z);
  });
  Field<double> expected(grid, [&](double x, double y, double z) {
    return std::exp(-3.0 * D * T) * std::sin(x) * std::sin(y) * std::sin(z);
  });

  diffuse(psi, D, dt, n_steps, RK4{});

  // RK4 has near-zero time error; the residual is dominated by the spatial
  // O(h^2) stencil error from Phase 1's Laplacian.
  const double rel_err = maxAbsDiff(psi, expected) / maxAbs(expected);
  EXPECT_LT(rel_err, 1e-2) << "relative max-norm error was " << rel_err;
}

TEST(IntegratorTest, CombinedRefinementOrder) {
  // Sweep N keeping physical T fixed. Both Euler and RK4 use dt at 40% of
  // their respective heat-equation CFL bounds. Both should show
  // log-log slope ~2 (combined error is spatial-dominated at CFL-limited dt).
  const double L = 2.0 * kPi;
  const double D = 0.1;
  const double T = 1.0;

  auto run_sweep = [&](auto tag) {
    using Tag = decltype(tag);
    const double cfl_limit = Tag::diffusionCFLLimit;
    const std::vector<int> Ns = {16, 32, 64};
    std::vector<double> errors;
    errors.reserve(Ns.size());
    for (int N : Ns) {
      const double h = L / static_cast<double>(N);
      const double inv_sum = 3.0 / (h * h);
      const double dt = 0.4 * cfl_limit / (D * inv_sum);
      const int n_steps = static_cast<int>(std::ceil(T / dt));
      const double T_actual = static_cast<double>(n_steps) * dt;

      Grid grid(N, N, N, Vec3d(h, h, h));
      Field<double> psi(grid, [&](double x, double y, double z) {
        return std::sin(x) * std::sin(y) * std::sin(z);
      });
      Field<double> expected(grid, [&](double x, double y, double z) {
        return std::exp(-3.0 * D * T_actual) * std::sin(x) * std::sin(y) * std::sin(z);
      });

      diffuse(psi, D, dt, n_steps, tag);
      errors.push_back(maxAbsDiff(psi, expected) / maxAbs(expected));
    }

    const std::size_t n = Ns.size();
    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      const double h = L / static_cast<double>(Ns[i]);
      const double x = std::log(h);
      const double y = std::log(errors[i]);
      sum_x += x; sum_y += y; sum_xy += x * y; sum_xx += x * x;
    }
    return (static_cast<double>(n) * sum_xy - sum_x * sum_y) /
           (static_cast<double>(n) * sum_xx - sum_x * sum_x);
  };

  const double order_euler = run_sweep(ExplicitEuler{});
  const double order_rk4 = run_sweep(RK4{});

  EXPECT_GT(order_euler, 1.8) << "Euler combined slope was " << order_euler;
  EXPECT_LT(order_euler, 2.5) << "Euler combined slope was " << order_euler;
  EXPECT_GT(order_rk4, 1.8) << "RK4 combined slope was " << order_rk4;
  EXPECT_LT(order_rk4, 2.5) << "RK4 combined slope was " << order_rk4;
}

TEST(IntegratorTest, RK4LessErrorThanEulerAtSameDt) {
  // Pure linear-decay ODE: dpsi/dt = -alpha * psi. No Laplacian, no spatial
  // discretization error -- so the entire error is from time integration.
  // For the diffusion problem, Euler's spatial-stencil error and time-
  // discretization error partially cancel at small dt; this test isolates
  // the time-discretization comparison.
  const int N = 4;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  const double alpha = 1.0;
  const double dt = 0.5;
  const int n_steps = 4;
  const double T = static_cast<double>(n_steps) * dt;

  Grid grid(N, N, N, Vec3d(h, h, h));
  auto make_psi = [&]() { return Field<double>(grid, 1.0); };
  auto rhs = [alpha](const Field<double>& f) {
    Field<double> out(f.getGrid());
    const std::size_t N_ = out.getSize();
    const double* p = f.data();
    double* o = out.data();
    for (std::size_t n = 0; n < N_; ++n) o[n] = -alpha * p[n];
    return out;
  };

  Field<double> psi_euler = make_psi();
  Field<double> psi_rk4 = make_psi();
  integrate(psi_euler, rhs, dt, n_steps, ExplicitEuler{});
  integrate(psi_rk4, rhs, dt, n_steps, RK4{});

  const double expected = std::exp(-alpha * T);
  const double err_euler = std::abs(psi_euler.data()[0] - expected);
  const double err_rk4 = std::abs(psi_rk4.data()[0] - expected);

  // RK4's per-step error is O((alpha*dt)^5/120) ~ 3e-4; Euler's is
  // O((alpha*dt)^2/2) ~ 1.25e-1. Ratio per step is ~400; integrated over
  // 4 steps the ratio remains around two orders of magnitude.
  EXPECT_LT(err_rk4 * 10.0, err_euler)
      << "err_euler=" << err_euler << " err_rk4=" << err_rk4;
}

TEST(IntegratorTest, RK4TimeAccuracyOrder4) {
  // Pure linear-decay ODE (no PDE, no spatial discretization) so the
  // entire error is RK4's time-discretization error. Roadmap acceptance
  // for Phase 6: "error scales as dt^4."
  //
  // RK4 on dpsi/dt = -alpha * psi yields the truncated Taylor series
  //   g(dt) = 1 - alpha*dt + (alpha*dt)^2/2 - (alpha*dt)^3/6 + (alpha*dt)^4/24
  // and after n = T/dt steps the global error is O(dt^4 * T).
  //
  // Sweep dt in {1.0, 0.5, 0.25} with alpha=1 and T=1. Halving dt should
  // reduce the error by ~16x. Log-log slope ~4.
  const int N = 8;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  const double alpha = 1.0;
  const double T = 1.0;

  Grid grid(N, N, N, Vec3d(h, h, h));
  auto make_psi = [&]() { return Field<double>(grid, 1.0); };
  auto rhs = [alpha](const Field<double>& f) {
    Field<double> out(f.getGrid());
    const std::size_t N_ = out.getSize();
    const double* p = f.data();
    double* o = out.data();
    for (std::size_t n = 0; n < N_; ++n) o[n] = -alpha * p[n];
    return out;
  };

  const double expected = std::exp(-alpha * T);
  const std::vector<double> dts = {1.0, 0.5, 0.25};
  std::vector<double> errors;
  errors.reserve(dts.size());

  for (double dt : dts) {
    const int n = static_cast<int>(std::round(T / dt));
    Field<double> psi = make_psi();
    integrate(psi, rhs, dt, n, RK4{});
    errors.push_back(std::abs(psi.data()[0] - expected));
  }

  // Halving dt should reduce error by ~16 (slope 4 globally).
  for (std::size_t i = 1; i < errors.size(); ++i) {
    const double ratio = errors[i - 1] / errors[i];
    EXPECT_GT(ratio, 8.0) << "ratio at step " << i << " was " << ratio
                          << " (errors: " << errors[i - 1] << " -> " << errors[i] << ")";
  }

  // Log-log slope.
  double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;
  for (std::size_t i = 0; i < dts.size(); ++i) {
    const double x = std::log(dts[i]);
    const double y = std::log(errors[i]);
    sum_x += x; sum_y += y; sum_xy += x * y; sum_xx += x * x;
  }
  const double n = static_cast<double>(dts.size());
  const double order = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);
  EXPECT_GT(order, 3.5) << "RK4 time-accuracy log-log slope was " << order;
  EXPECT_LT(order, 5.0) << "RK4 time-accuracy log-log slope was " << order;
}

TEST(IntegratorTest, RK4StabilityGap) {
  // Pick D, dt, h such that the CFL number lies between Euler's 0.5 and
  // RK4's 0.6963 -- Euler should throw, RK4 should run successfully.
  const int N = 16;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  const double D = 0.1;
  const double inv_sum = 3.0 / (h * h);
  // Target CFL = 0.6 (between 0.5 and 0.6963).
  const double dt = 0.6 / (D * inv_sum);
  const double cfl = D * dt * inv_sum;
  ASSERT_GT(cfl, 0.5);
  ASSERT_LT(cfl, 0.6963);

  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi_euler(grid, 1.0);
  Field<double> psi_rk4(grid, 1.0);

  EXPECT_THROW(diffuse(psi_euler, D, dt, 1, ExplicitEuler{}), std::invalid_argument);
  EXPECT_NO_THROW(diffuse(psi_rk4, D, dt, 1, RK4{}));
}

TEST(IntegratorTest, GenericIntegrateDispatchEulerEqualsExplicitEuler) {
  const int N = 16;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  const double dt = 0.001;
  const int n_steps = 50;

  Grid grid(N, N, N, Vec3d(h, h, h));
  auto make_psi = [&]() {
    return Field<double>(grid, [&](double x, double y, double z) {
      return std::sin(x) * std::cos(y) * std::sin(z);
    });
  };

  auto rhs = [](const Field<double>& f) {
    Field<double> out(f.getGrid());
    const std::size_t N_ = out.getSize();
    const double* p = f.data();
    double* o = out.data();
    for (std::size_t n = 0; n < N_; ++n) o[n] = -0.1 * p[n];
    return out;
  };

  Field<double> psi_a = make_psi();
  Field<double> psi_b = make_psi();
  explicitEuler(psi_a, rhs, dt, n_steps);
  integrate(psi_b, rhs, dt, n_steps, ExplicitEuler{});

  EXPECT_TRUE(bitEqualField(psi_a, psi_b))
      << "explicitEuler() and integrate(... ExplicitEuler{}) diverged";
}
