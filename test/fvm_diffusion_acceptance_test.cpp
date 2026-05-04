// Phase 14 — heterogeneous-D diffusion end-to-end acceptance test.
//
// Discharges the roadmap acceptance bullet
//   "Heterogeneous-D diffusion of a known closed-form decay matches the
//    analytical answer at t = 1 to O(h^2); plus correct qualitative
//    behavior on genuinely heterogeneous D."
//
// **Why two tests, not one.** A clean closed-form solution to the
// continuous PDE $\partial_t\psi = \nabla\cdot(D(\mathbf x)\nabla\psi)$ is
// only available in special cases. The naive manufactured solution
// "psi = exp(-D(x) t) sin(z)" does not actually satisfy the PDE for
// t > 0 — once psi acquires x-dependence from the position-dependent
// decay rate, the x-axis flux $\partial_x[D \partial_x \psi]$ becomes
// non-zero and the analytical reference drifts away from the discrete
// trajectory. The two tests below combine to discharge the acceptance
// bar:
//
//   1. `RecoversAnalyticalEigenfunctionAtOrder2` — runs the
//      heterogeneous-D solve::diffuse code path with a uniform D
//      field. The trajectory is the well-known
//      psi(t) = exp(-3 D t) sin(x) sin(y) sin(z), and convergence is
//      order 2 in h. This exercises every line of the new path
//      (`cellLaplacian` with harmonic-mean averaging, the
//      heterogeneous-D `solve::diffuse` overload, the CFL check, the
//      D-positivity validation, the integrator loop) on a problem
//      with a known closed-form answer.
//
//   2. `HeterogeneousDQualitativeProperties` — runs with a genuinely
//      varying D(x) and a Gaussian IC; asserts the three properties
//      that hold for ANY D > 0: total mass is conserved (FVM
//      conservation + periodic BCs); max|psi| is monotonically
//      non-increasing across snapshots (max-principle / dissipativity
//      under symmetric heterogeneous diffusion); psi >= 0 throughout
//      (positivity preservation for a positive IC).

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/func/Integrate.hpp>
#include <gridcalc/solve/Diffusion.hpp>

namespace {

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
using gridcalc::func::integrate;

constexpr double kPi = 3.14159265358979323846;

// =============================================================================
// 1. Order-2 convergence on the heterogeneous-D path with uniform D.
// =============================================================================

double eigenfunctionErrorAtN(int N, double D_value, double T) {
    const double L = 2.0 * kPi;
    const double h = L / static_cast<double>(N);
    Grid g(N, N, N, Vec3d(h, h, h));

    Field<double> psi(g, [](double x, double y, double z) {
        return std::sin(x) * std::sin(y) * std::sin(z);
    });
    Field<double> D(g, D_value);

    // RK4 CFL: D * dt * (3 / h^2) <= 0.6963.
    const double dt_limit = 0.6963 * h * h / (3.0 * D_value);
    const double dt = 0.4 * dt_limit;
    const int n_steps = static_cast<int>(std::ceil(T / dt));
    const double t_actual = static_cast<double>(n_steps) * dt;

    gridcalc::solve::diffuse(psi, D, dt, n_steps, gridcalc::solve::RK4{});

    // Analytical reference: psi(t) = exp(-3 D t) sin(x) sin(y) sin(z).
    const double decay = std::exp(-3.0 * D_value * t_actual);
    double err = 0.0;
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                const double x = static_cast<double>(i) * h;
                const double y = static_cast<double>(j) * h;
                const double z = static_cast<double>(k) * h;
                const double psi_ref = decay * std::sin(x) * std::sin(y) * std::sin(z);
                err = std::max(err, std::abs(psi(i, j, k) - psi_ref));
            }
        }
    }
    return err;
}

TEST(FvmDiffusionAcceptanceTest, RecoversAnalyticalEigenfunctionAtOrder2) {
    const std::vector<int> Ns = {16, 24, 32};
    const double D_value = 0.1;
    const double T = 1.0;

    std::vector<double> errors;
    errors.reserve(Ns.size());
    for (int N : Ns) {
        errors.push_back(eigenfunctionErrorAtN(N, D_value, T));
    }

    // Slopes between consecutive (h, error) pairs. Expect ~2 within a
    // generous [1.6, 2.4] band.
    for (std::size_t i = 0; i + 1 < Ns.size(); ++i) {
        const double slope = std::log2(errors[i] / errors[i + 1]) /
                             std::log2(static_cast<double>(Ns[i + 1]) /
                                       static_cast<double>(Ns[i]));
        EXPECT_GT(slope, 1.6) << "N=" << Ns[i] << "->" << Ns[i + 1]
                              << ": err=" << errors[i] << " -> " << errors[i + 1]
                              << ", slope=" << slope;
        EXPECT_LT(slope, 2.4) << "N=" << Ns[i] << "->" << Ns[i + 1]
                              << ": err=" << errors[i] << " -> " << errors[i + 1]
                              << ", slope=" << slope;
    }
}

// =============================================================================
// 2. Qualitative properties on a genuinely heterogeneous D field.
// =============================================================================

TEST(FvmDiffusionAcceptanceTest, HeterogeneousDQualitativeProperties) {
    constexpr int N = 24;
    const double L = 2.0 * kPi;
    const double h = L / N;
    Grid g(N, N, N, Vec3d(h, h, h));

    // Gaussian peak IC, centered at (pi, pi, pi).
    const double sigma = 0.6;
    Field<double> psi(g, [&](double x, double y, double z) {
        const double dx = x - kPi;
        const double dy = y - kPi;
        const double dz = z - kPi;
        return std::exp(-(dx * dx + dy * dy + dz * dz) / (2.0 * sigma * sigma));
    });
    // D varies in x: ranges over [0.05, 0.15], strictly positive.
    Field<double> D(g, [](double x, double, double) {
        return 0.10 + 0.05 * std::cos(x);
    });

    const double mass_initial = integrate(psi);
    auto maxAbs = [](const Field<double>& f) {
        double m = 0.0;
        for (double v : f) m = std::max(m, std::abs(v));
        return m;
    };
    const double peak_initial = maxAbs(psi);

    // RK4 step with D_max = 0.15.
    const double D_max = 0.15;
    const double dt = 0.4 * 0.6963 * h * h / (3.0 * D_max);
    const int n_steps_per_chunk = 50;
    const int n_chunks = 4;

    double prev_peak = peak_initial;
    double min_psi = 0.0;
    for (int chunk = 0; chunk < n_chunks; ++chunk) {
        gridcalc::solve::diffuse(psi, D, dt, n_steps_per_chunk, gridcalc::solve::RK4{});

        // Mass conservation across every chunk.
        const double mass_now = integrate(psi);
        EXPECT_NEAR(mass_now, mass_initial, 1e-10 * std::abs(mass_initial))
            << "mass drift after chunk " << chunk;

        // max|psi| monotonically non-increasing (dissipativity).
        const double peak_now = maxAbs(psi);
        EXPECT_LE(peak_now, prev_peak + 1e-12)
            << "max|psi| grew between chunks " << chunk - 1 << " and " << chunk
            << ": " << prev_peak << " -> " << peak_now;
        prev_peak = peak_now;

        // psi stays non-negative everywhere.
        for (double v : psi) min_psi = std::min(min_psi, v);
    }
    EXPECT_GE(min_psi, -1e-12) << "psi went negative; observed min = " << min_psi;
}

}  // namespace
