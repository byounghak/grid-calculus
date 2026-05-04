#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/func/Integrate.hpp>
#include <gridcalc/solve/Diffusion.hpp>
#include <gridcalc/solve/ExplicitEuler.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
using gridcalc::func::integrate;
using gridcalc::solve::diffuse;
using gridcalc::solve::explicitEuler;

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

}  // namespace

TEST(DiffusionTest, TrigEigenfunctionDecaysToAnalytic) {
  // psi_0 = sin(x) sin(y) sin(z) is an eigenfunction of the Laplacian on
  // [0, 2*pi]^3 with eigenvalue -3, so psi(t) = exp(-3 D t) * psi_0.
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

  diffuse(psi, D, dt, n_steps);

  // Discrete eigenvalue is -2(1 - cos(h))/h^2 per axis, summed over 3 axes.
  // Total leading-order error scales as h^2; at N=32 this is a few percent.
  const double rel_err = maxAbsDiff(psi, expected) / maxAbs(expected);
  EXPECT_LT(rel_err, 5e-2) << "relative max-norm error was " << rel_err;
}

TEST(DiffusionTest, TrigEigenfunctionSecondOrderConvergence) {
  // Sweep N keeping the physical end time T fixed. Because the explicit-Euler
  // CFL forces dt ~ h^2, the time-discretization error inherits the same
  // O(h^2) scaling as the spatial error -- so the combined error is O(h^2).
  const double L = 2.0 * kPi;
  const double D = 0.1;
  const double T = 1.0;

  const std::vector<int> Ns = {16, 32, 64};
  std::vector<double> errors;
  errors.reserve(Ns.size());

  for (int N : Ns) {
    const double h = L / static_cast<double>(N);
    // Stay well under the CFL limit dt = h^2 / (6 D).
    const double dt_limit = h * h / (6.0 * D);
    const double dt = 0.4 * dt_limit;
    const int n_steps = static_cast<int>(std::ceil(T / dt));
    const double T_actual = static_cast<double>(n_steps) * dt;

    Grid grid(N, N, N, Vec3d(h, h, h));
    Field<double> psi(grid, [&](double x, double y, double z) {
      return std::sin(x) * std::sin(y) * std::sin(z);
    });
    Field<double> expected(grid, [&](double x, double y, double z) {
      return std::exp(-3.0 * D * T_actual) * std::sin(x) * std::sin(y) * std::sin(z);
    });

    diffuse(psi, D, dt, n_steps);
    errors.push_back(maxAbsDiff(psi, expected) / maxAbs(expected));
  }

  // Log-log slope of error vs h. Allow a generous upper bound of 2.5 because
  // explicit-Euler's O(dt) time error rides on top of the spatial O(h^2)
  // and dt ~ h^2 makes the time error look superlinear in h.
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
  EXPECT_LT(order, 2.5) << "log-log convergence order was " << order;
}

TEST(DiffusionTest, GaussianSanity) {
  // Gaussian IC; assert finite values, mass conservation, and peak decay.
  const int N = 64;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  const double D = 0.05;
  const double dt = 0.4 * h * h / (6.0 * D);
  const int n_steps = 100;

  Grid grid(N, N, N, Vec3d(h, h, h));
  const double sigma = 0.5;
  Field<double> psi(grid, [&](double x, double y, double z) {
    const double dx = x - kPi;
    const double dy = y - kPi;
    const double dz = z - kPi;
    return std::exp(-(dx * dx + dy * dy + dz * dz) / (2.0 * sigma * sigma));
  });

  const double mass_before = integrate(psi);
  const double peak_before = maxAbs(psi);

  diffuse(psi, D, dt, n_steps);

  // No NaN / Inf.
  for (double v : psi) {
    EXPECT_TRUE(std::isfinite(v));
  }

  // Mass is conserved on a periodic domain (boundary flux is zero).
  const double mass_after = integrate(psi);
  EXPECT_NEAR(mass_after, mass_before, 1e-10 * std::abs(mass_before))
      << "mass_before=" << mass_before << " mass_after=" << mass_after;

  // Gaussian has spread out -- peak height is strictly smaller.
  const double peak_after = maxAbs(psi);
  EXPECT_LT(peak_after, peak_before)
      << "peak_before=" << peak_before << " peak_after=" << peak_after;
}

TEST(DiffusionTest, CFLViolationThrows) {
  const int N = 32;
  const double h = 2.0 * kPi / static_cast<double>(N);
  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, 1.0);

  // CFL bound: D * dt * (3 / h^2) <= 0.5  ->  D * dt <= h^2 / 6.
  // Pick D * dt that overshoots by 2x.
  const double D = 1.0;
  const double dt = h * h / 3.0;  // -> CFL number = 1.0, exceeds 0.5.

  EXPECT_THROW(diffuse(psi, D, dt, 1), std::invalid_argument);
}

TEST(DiffusionTest, GenericExplicitEulerOnZeroRhs) {
  const int N = 8;
  const double h = 2.0 * kPi / static_cast<double>(N);
  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [&](double x, double y, double z) {
    return std::sin(x) * std::sin(y) * std::sin(z);
  });

  // Snapshot psi values to verify bit-identity after zero-RHS integration.
  std::vector<double> snapshot(psi.data(), psi.data() + psi.getSize());

  auto zero_rhs = [](const Field<double>& f) {
    return Field<double>(f.getGrid(), 0.0);
  };

  explicitEuler(psi, zero_rhs, 0.1, 100);

  for (std::size_t n = 0; n < psi.getSize(); ++n) {
    EXPECT_EQ(psi.data()[n], snapshot[n]);
  }
}

TEST(DiffusionTest, RejectsNegativeNSteps) {
  const int N = 4;
  const double h = 0.1;
  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, 0.0);

  EXPECT_THROW(diffuse(psi, 0.1, 0.001, -1), std::invalid_argument);

  auto rhs = [](const Field<double>& f) { return Field<double>(f.getGrid(), 0.0); };
  EXPECT_THROW(explicitEuler(psi, rhs, 0.001, -1), std::invalid_argument);
}

// =============================================================================
// Heterogeneous-D solve::diffuse overload (Phase 14, FVM).
// =============================================================================

TEST(DiffusionTest, HeterogeneousD_AgreesWithConstantDOnUniformField) {
  // With a uniform D-field, the heterogeneous-D path must produce the
  // same trajectory as the constant-D path to round-off (FVM reduces to
  // FD on uniform grids with uniform D, modulo summation order).
  const int N = 16;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  const double D_const = 0.1;
  const double dt = 0.4 * h * h / (6.0 * D_const);
  const int n_steps = 50;

  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi_const(grid, [&](double x, double y, double z) {
    return std::sin(x) * std::cos(y) * std::sin(z);
  });
  Field<double> psi_het = psi_const;  // identical IC
  Field<double> D_field(grid, D_const);

  diffuse(psi_const, D_const, dt, n_steps);
  diffuse(psi_het, D_field, dt, n_steps);

  EXPECT_LT(maxAbsDiff(psi_const, psi_het), 1e-12);
}

TEST(DiffusionTest, HeterogeneousD_RejectsCFLViolation) {
  const int N = 32;
  const double h = 2.0 * kPi / static_cast<double>(N);
  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, 1.0);
  // D-field with high spots: D_max = 1.0; CFL needs D_max * dt * 3/h^2 <= 0.5
  // -> dt <= h^2/6. Overshoot.
  Field<double> D(grid, [](double x, double, double) { return 0.5 + 0.5 * std::cos(x); });
  // smallest positive D at cos = -1 is 0.0 — bump to ensure D > 0:
  for (std::size_t idx = 0; idx < D.getSize(); ++idx) {
    if (D.data()[idx] <= 0.0) D.data()[idx] = 1e-6;
  }
  const double dt = h * h / 3.0;  // CFL = 1.0 with D_max ~ 1.0, exceeds 0.5.

  EXPECT_THROW(diffuse(psi, D, dt, 1), std::invalid_argument);
}

TEST(DiffusionTest, HeterogeneousD_RejectsNonpositiveD) {
  const int N = 8;
  const double h = 0.1;
  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, 1.0);

  Field<double> D(grid, 1.0);
  D(3, 4, 5) = -0.01;  // single bad cell

  EXPECT_THROW(diffuse(psi, D, 1e-4, 1), std::invalid_argument);

  D(3, 4, 5) = 0.0;  // boundary case: zero is also rejected
  EXPECT_THROW(diffuse(psi, D, 1e-4, 1), std::invalid_argument);
}

TEST(DiffusionTest, HeterogeneousD_MassConservedOverTimeStepping) {
  // Mass conservation: the spatial sum of psi must be constant across
  // many integration steps (FVM's per-step conservation property carries
  // through the integrator).
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  Grid grid(N, N, N, Vec3d(h, h, h));

  // Gaussian peak IC.
  const double sigma = 0.5;
  Field<double> psi(grid, [&](double x, double y, double z) {
    const double dx = x - kPi;
    const double dy = y - kPi;
    const double dz = z - kPi;
    return std::exp(-(dx * dx + dy * dy + dz * dz) / (2.0 * sigma * sigma));
  });
  Field<double> D(grid, [](double x, double, double) {
    return 0.05 + 0.04 * std::sin(x);  // ranges over [0.01, 0.09]
  });

  const double mass_before = integrate(psi);
  // CFL: D_max ~ 0.09; dt <= h^2 / (6 * D_max) for ExplicitEuler.
  const double D_max = 0.09;
  const double dt = 0.4 * h * h / (6.0 * D_max);
  const int n_steps = 100;
  diffuse(psi, D, dt, n_steps);
  const double mass_after = integrate(psi);

  EXPECT_NEAR(mass_after, mass_before, 1e-10 * std::abs(mass_before));
}
