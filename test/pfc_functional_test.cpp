#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/diff/Biharmonic.hpp>
#include <gridcalc/diff/Laplacian.hpp>
#include <gridcalc/func/Functional.hpp>
#include <gridcalc/func/Integrate.hpp>

#ifdef GRIDCALC_USE_FFT
#include <gridcalc/spectral/Derivatives.hpp>
#endif

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
using gridcalc::func::evaluate;
using gridcalc::func::integrate;

namespace {

constexpr double kPi = 3.14159265358979323846;

double piCubed() {
  return kPi * kPi * kPi;
}

// PFC integrand for the Phase 11 worked example. Used by both the production
// path (passed as a 4-arg lambda to func::evaluate) and by the spectral
// cross-check (called point-by-point with spectrally computed derivatives).
//
// f(psi, lap, biharm) = (1/2)(q0^4 - eps) * psi^2
//                     + q0^2 * psi * lap
//                     + (1/2) * psi * biharm
//                     + (1/4) * psi^4
//
// which equals (1/2) psi [(q0^2 + nabla^2)^2 - eps] psi + (1/4) psi^4
// after expansion.
double pfcDensity(double psi, double lap, double biharm, double q0, double eps) {
  const double q0_sq = q0 * q0;
  const double q0_4 = q0_sq * q0_sq;
  const double psi_sq = psi * psi;
  return 0.5 * (q0_4 - eps) * psi_sq + q0_sq * psi * lap + 0.5 * psi * biharm +
         0.25 * psi_sq * psi_sq;
}

// Hand-computed reference for psi = sin(x) sin(y) sin(z) on [0, 2pi]^3:
//   integral psi^2     = pi^3
//   integral psi^4     = (3/4)^3 * pi^3 = 27 pi^3 / 64
//   lap psi            = -3 psi   =>  integral psi * lap     = -3 pi^3
//   biharm psi         =  9 psi   =>  integral psi * biharm  =  9 pi^3
// With q0 = 1, eps = 0.1:
//   F_ref = [(1/2)(1 - 0.1) + (-3) + (1/2)(9) + 27/256] * pi^3
//         = (0.45 - 3 + 4.5 + 0.10546875) * pi^3
//         = 2.05546875 * pi^3
double pfcReferenceValue(double q0, double eps) {
  const double q0_sq = q0 * q0;
  const double q0_4 = q0_sq * q0_sq;
  const double term_psi2 = 0.5 * (q0_4 - eps);
  const double term_lap = q0_sq * (-3.0);
  const double term_biharm = 0.5 * 9.0;
  const double term_psi4 = (1.0 / 4.0) * (27.0 / 64.0);
  return (term_psi2 + term_lap + term_biharm + term_psi4) * piCubed();
}

}  // namespace

TEST(PfcFunctionalTest, HandComputedAtN64) {
  // Roadmap acceptance bar: PFC functional value matches a hand-computed
  // reference for psi = sin(x)sin(y)sin(z) on [0, 2pi]^3 within stencil-order
  // error. At N=64, the dominant discrete error comes from the biharmonic
  // term (Phase 8 fused stencil at default Order=2, leading term O(h^2)).
  const int N = 64;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  const double q0 = 1.0;
  const double eps = 0.1;

  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [](double x, double y, double z) {
    return std::sin(x) * std::sin(y) * std::sin(z);
  });

  auto pfc_lambda = [q0, eps](double p, const Vec3d&, double lp, double bp) {
    return pfcDensity(p, lp, bp, q0, eps);
  };
  const double F = evaluate(psi, pfc_lambda);
  const double F_ref = pfcReferenceValue(q0, eps);

  const double rel_err = std::abs(F - F_ref) / std::abs(F_ref);
  EXPECT_LT(rel_err, 2e-2) << "F=" << F << " F_ref=" << F_ref << " rel_err=" << rel_err;
}

TEST(PfcFunctionalTest, SecondOrderConvergence) {
  // Confirms the FD path (4-arg arity, default Order=2 derivatives) converges
  // at the documented O(h^2) rate over N in {16, 32, 64}.
  const double L = 2.0 * kPi;
  const double q0 = 1.0;
  const double eps = 0.1;
  const double F_ref = pfcReferenceValue(q0, eps);

  const std::vector<int> Ns = {16, 32, 64};
  std::vector<double> errors;
  errors.reserve(Ns.size());

  for (int N : Ns) {
    const double h = L / static_cast<double>(N);
    Grid grid(N, N, N, Vec3d(h, h, h));
    Field<double> psi(grid, [](double x, double y, double z) {
      return std::sin(x) * std::sin(y) * std::sin(z);
    });
    auto pfc_lambda = [q0, eps](double p, const Vec3d&, double lp, double bp) {
      return pfcDensity(p, lp, bp, q0, eps);
    };
    const double F = evaluate(psi, pfc_lambda);
    errors.push_back(std::abs(F - F_ref) / std::abs(F_ref));
  }

  for (std::size_t i = 1; i < errors.size(); ++i) {
    const double ratio = errors[i - 1] / errors[i];
    EXPECT_GT(ratio, 3.5) << "halving ratio at step " << i << " was " << ratio;
    EXPECT_LT(ratio, 4.5) << "halving ratio at step " << i << " was " << ratio;
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
  EXPECT_GT(order, 1.7) << "log-log convergence order was " << order;
  EXPECT_LT(order, 2.3) << "log-log convergence order was " << order;
}

TEST(PfcFunctionalTest, PsiOnlyArityStillWorks) {
  // Regression check: the new 4-arg branch must not shadow the 1-arg branch.
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);

  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [](double x, double y, double z) {
    return std::sin(x) * std::sin(y) * std::sin(z);
  });

  auto f = [](double p) { return p * p; };
  const double F = evaluate(psi, f);
  EXPECT_NEAR(F, piCubed(), 1e-12);
}

TEST(PfcFunctionalTest, GinzburgLandauStillWorks) {
  // Regression check: the existing Phase 4 GL hand-computed value reproduces
  // via the 3-arg dispatch path (squeezed between the new 4-arg branch above
  // and the 2-arg branch below in the if-constexpr ladder).
  const int N = 64;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  const double W = 1.0;

  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [](double x, double y, double z) {
    return std::sin(x) * std::sin(y) * std::sin(z);
  });

  auto gl_density = [W](double p, const Vec3d& g, double) {
    const double dpsi2 = p * p - 1.0;
    return 0.5 * g.squaredNorm() + W * dpsi2 * dpsi2;
  };

  const double F = evaluate(psi, gl_density);
  const double F_ref = (3.0 / 2.0 + 411.0 / 64.0) * piCubed();

  const double rel_err = std::abs(F - F_ref) / std::abs(F_ref);
  EXPECT_LT(rel_err, 1e-2) << "F=" << F << " F_ref=" << F_ref;
}

#ifdef GRIDCALC_USE_FFT

TEST(PfcFunctionalTest, SpectralCrossCheck) {
  // Independent numerical path: materialize lap and biharm via the Phase 9
  // spectral operators (exact for band-limited inputs), assemble the PFC
  // integrand by hand, reduce with func::integrate. Since psi is a single
  // Fourier mode, the spectral derivatives are exact to machine precision and
  // the resulting F must match the analytical reference to ~1e-12 relative.
  // This confirms the analytical reference itself, independently of the FD
  // biharmonic stencil exercised by HandComputedAtN64.
  const int N = 32;
  const double L = 2.0 * kPi;
  const double h = L / static_cast<double>(N);
  const double q0 = 1.0;
  const double eps = 0.1;

  Grid grid(N, N, N, Vec3d(h, h, h));
  Field<double> psi(grid, [](double x, double y, double z) {
    return std::sin(x) * std::sin(y) * std::sin(z);
  });

  const auto lap = gridcalc::spectral::laplacian(psi);
  const auto biharm = gridcalc::spectral::biharmonic(psi);

  Field<double> integrand(grid);
  for (int k = 0; k < N; ++k) {
    for (int j = 0; j < N; ++j) {
      for (int i = 0; i < N; ++i) {
        integrand(i, j, k) =
            pfcDensity(psi(i, j, k), lap(i, j, k), biharm(i, j, k), q0, eps);
      }
    }
  }
  const double F_spectral = integrate(integrand);
  const double F_ref = pfcReferenceValue(q0, eps);

  const double rel_err = std::abs(F_spectral - F_ref) / std::abs(F_ref);
  EXPECT_LT(rel_err, 1e-12)
      << "F_spectral=" << F_spectral << " F_ref=" << F_ref << " rel_err=" << rel_err;
}

#endif  // GRIDCALC_USE_FFT
