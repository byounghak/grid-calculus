// Acceptance + smoke tests for the Cahn-Hilliard demo.
//
// Roadmap acceptance bar (specs/roadmap.md, Phase 12): "Coarsening
// dynamics qualitatively match published CH benchmarks." We discharge it
// in CI with three programmatic checks on a small simulation:
//
//   * `LyapunovDecay` — total free energy F(t) is monotonically
//     non-increasing across snapshots (the dF/dt = -M ∫|grad mu|^2 ≤ 0
//     property of CH dynamics).
//
//   * `MonotoneCoarsening` — gradient energy ∫|grad psi|^2 dV decreases
//     between a pre-coarsening snapshot (at the gradient-energy peak,
//     just after interfaces saturate) and a post-coarsening snapshot
//     (one decade further in time). On the chosen parameters the drop
//     is ~5%–8% across multiple seeds; the test asserts a 3% drop.
//
//   * `MassConservation` — the spatial mean is conserved at
//     round-off across the run.
//
// Plus three quick correctness tests on the RHS itself
// (`ZeroRhsForConstantField`, `RhsIsConservative`,
// `SeedReproducibility`).
//
// Total wall-clock target: < 60 s on Ubuntu Clang debug. Local
// measurement: ~30 s on Apple Clang (clang-debug, -O0 -g) for the
// 2500-step coarsening run.

#include <common/cahn_hilliard.hpp>

#include <cmath>
#include <cstdint>

#include <gtest/gtest.h>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/solve/Integrator.hpp>

namespace {

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
namespace ch = gridcalc::examples::cahn_hilliard;

TEST(CahnHilliardTest, ZeroRhsForConstantField) {
    Grid g(8, 8, 8, Vec3d(1.0, 1.0, 1.0));
    Field<double> psi(g, 0.42);
    auto rhs = ch::computeRhs(psi, /*M=*/1.0, /*kappa=*/1.0);

    const std::size_t n = rhs.getSize();
    double max_abs = 0.0;
    for (std::size_t idx = 0; idx < n; ++idx) {
        max_abs = std::max(max_abs, std::abs(rhs.data()[idx]));
    }
    EXPECT_LT(max_abs, 1e-12);
}

TEST(CahnHilliardTest, RhsIsConservative) {
    // d_t psi = M lap(mu); on a periodic domain, sum of lap(mu) over all
    // grid points is zero to round-off (telescoping central differences).
    // Equivalent to total mass being conserved by the dynamics.
    Grid g(16, 16, 16, Vec3d(1.0, 1.0, 1.0));
    Field<double> psi = ch::makeRandomIc(g, /*seed=*/123, /*amplitude=*/0.3);

    auto rhs = ch::computeRhs(psi, /*M=*/1.0, /*kappa=*/1.0);
    double sum = 0.0;
    const std::size_t n = rhs.getSize();
    for (std::size_t idx = 0; idx < n; ++idx) sum += rhs.data()[idx];

    // Reference scale for the absolute tolerance: max |rhs| × N points.
    double max_abs = 0.0;
    for (std::size_t idx = 0; idx < n; ++idx) {
        max_abs = std::max(max_abs, std::abs(rhs.data()[idx]));
    }
    EXPECT_LT(std::abs(sum), 1e-10 * max_abs * static_cast<double>(n));
}

TEST(CahnHilliardTest, SeedReproducibility) {
    Grid g(8, 8, 8, Vec3d(1.0, 1.0, 1.0));
    auto a = ch::makeRandomIc(g, /*seed=*/777, /*amplitude=*/0.1);
    auto b = ch::makeRandomIc(g, /*seed=*/777, /*amplitude=*/0.1);
    ASSERT_EQ(a.getSize(), b.getSize());
    for (std::size_t idx = 0; idx < a.getSize(); ++idx) {
        EXPECT_EQ(a.data()[idx], b.data()[idx]) << "mismatch at flat index " << idx;
    }
}

// One simulation produces all three coarsening-property assertions —
// the run itself is the expensive part. Test parameters are chosen so
// that:
//   * the linear (spinodal) regime ends by t ~ 6 (initial amplitude
//     0.1, sigma_max = M/(4*kappa) = 0.5),
//   * gradient energy E_grad peaks near t ~ 30 (interfaces saturate to
//     |psi|~1 and the system enters the coarsening regime),
//   * by t = 50 (1000 RK4 steps later) E_grad has dropped 5%–8% across
//     at least seeds {1, 2, 7, 42, 1234} (locally measured); the test
//     asserts a 3% drop, which is conservative against
//     uniform_real_distribution differences across libstdc++/libc++.
class CahnHilliardCoarseningRun : public ::testing::Test {
 protected:
    static constexpr int kNx = 24;
    static constexpr double kDt = 0.02;
    static constexpr double kKappa = 0.5;
    static constexpr double kMobility = 1.0;
    static constexpr double kAmplitude = 0.1;
    static constexpr std::uint64_t kSeed = 42;
    // Snapshots at t = 0, 10, 20, 30, 50.
    static constexpr int kStepsPerSegment[] = {500, 500, 500, 1000};

    struct Snapshot {
        int step;
        double t;
        double F;
        double E_grad;
        double mean;
    };

    // Run once per process and cache; the three TEST_Fs in this fixture
    // share the result so the ~30 s simulation does not repeat per test.
    static const std::vector<Snapshot>& Snaps() {
        static const std::vector<Snapshot> snaps = []() {
            Grid grid(kNx, kNx, kNx, Vec3d(1.0, 1.0, 1.0));
            Field<double> psi = ch::makeRandomIc(grid, kSeed, kAmplitude);
            auto rhs = [](const Field<double>& y) { return ch::computeRhs(y, kMobility, kKappa); };

            std::vector<Snapshot> out;
            out.push_back({0, 0.0, ch::computeFreeEnergy(psi, kKappa),
                           ch::computeGradientEnergy(psi), ch::computeMean(psi)});

            int step = 0;
            for (int seg = 0; seg < 4; ++seg) {
                const int seg_steps = kStepsPerSegment[seg];
                gridcalc::solve::integrate(psi, rhs, kDt, seg_steps, gridcalc::solve::RK4{});
                step += seg_steps;
                out.push_back({step, static_cast<double>(step) * kDt,
                               ch::computeFreeEnergy(psi, kKappa),
                               ch::computeGradientEnergy(psi), ch::computeMean(psi)});
            }
            return out;
        }();
        return snaps;
    }
};

constexpr int CahnHilliardCoarseningRun::kStepsPerSegment[];

TEST_F(CahnHilliardCoarseningRun, LyapunovDecay) {
    const auto& snaps = Snaps();
    ASSERT_EQ(snaps.size(), 5u);
    // F is the canonical CH Lyapunov functional. Allow a tiny
    // tolerance for explicit-RK4 round-off; the actual drop per
    // segment is hundreds, so 1e-6 of |F[0]| is far below signal.
    const double tol_F = 1e-6 * std::abs(snaps[0].F) + 1e-9;
    for (std::size_t i = 0; i + 1 < snaps.size(); ++i) {
        EXPECT_LE(snaps[i + 1].F, snaps[i].F + tol_F)
            << "F should be non-increasing; snap " << i << " -> " << (i + 1)
            << ": " << snaps[i].F << " -> " << snaps[i + 1].F;
    }
}

TEST_F(CahnHilliardCoarseningRun, MonotoneCoarsening) {
    const auto& snaps = Snaps();
    // snap[3] is at t=30 (gradient-energy peak); snap[4] is at t=50
    // (one decade into coarsening). The post-snapshot must be at least
    // 3% below the pre-snapshot — a robust margin given the locally
    // observed 5%-8% drop across seeds {1, 2, 7, 42, 1234}.
    const double E_pre = snaps[3].E_grad;
    const double E_post = snaps[4].E_grad;
    const double drop_fraction = (E_pre - E_post) / E_pre;
    EXPECT_GT(drop_fraction, 0.03)
        << "Gradient energy should decrease by at least 3% from t=30 (peak) "
           "to t=50 (coarsening). Observed E_pre="
        << E_pre << ", E_post=" << E_post << ", drop=" << (drop_fraction * 100.0) << "%.";
}

TEST_F(CahnHilliardCoarseningRun, MassConservation) {
    const auto& snaps = Snaps();
    const double mean_t0 = snaps.front().mean;
    const double mean_tf = snaps.back().mean;
    EXPECT_LT(std::abs(mean_tf - mean_t0), 1e-10)
        << "<psi> drifted: " << mean_t0 << " -> " << mean_tf;
}

}  // namespace
