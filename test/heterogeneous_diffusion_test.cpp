// CI acceptance gate for the heterogeneous-D diffusion demo.
//
// Runs a stripped-down version of the same simulation the demo binary
// runs (Gaussian IC + smooth heterogeneous D + RK4 explicit) once per
// process via a static-cached lambda; three TEST_F bodies share the
// result so the simulation does not repeat per test. Mirrors the
// Phase 12 `test/cahn_hilliard_test.cpp` pattern.
//
// Properties asserted (all hold for any `D > 0`):
//   1. Mass conservation: `∫ψ dV` is constant to round-off across
//      every snapshot.
//   2. L²-norm decay: `∫ψ² dV` is monotonically non-increasing
//      between consecutive snapshots (the textbook entropy-like
//      Lyapunov property `d/dt ∫ψ²/2 = -∫D|∇ψ|² ≤ 0`).
//   3. Positivity preservation: `ψ ≥ 0` at every grid point at every
//      snapshot (a positive Gaussian IC stays positive under
//      diffusion with D > 0).

#include <common/heterogeneous_diffusion.hpp>

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/solve/Diffusion.hpp>

namespace {

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
namespace hd = gridcalc::examples::heterogeneous_diffusion;

constexpr double kPi = 3.14159265358979323846;

class HeterogeneousDiffusionRun : public ::testing::Test {
 protected:
    static constexpr int kNx = 24;
    static constexpr double kD0 = 0.10;
    static constexpr double kEps = 0.05;
    static constexpr double kSigma = 0.6;
    static constexpr int kStepsPerSegment = 50;
    static constexpr int kNumSegments = 4;

    struct Snapshot {
        int step;
        double mass;
        double l2sq;
        double psi_min;
        double psi_max;
    };

    static const std::vector<Snapshot>& Snaps() {
        static const std::vector<Snapshot> snaps = []() {
            const double L = 2.0 * kPi;
            const double h = L / kNx;
            Grid grid(kNx, kNx, kNx, Vec3d(h, h, h));
            const Vec3d center(L / 2.0, L / 2.0, L / 2.0);

            Field<double> psi = hd::makeGaussianIc(grid, kSigma, center);
            Field<double> D = hd::makeHeterogeneousD(grid, kD0, kEps, 'x');

            // Auto-pick dt from the RK4 CFL bound.
            const double D_max = kD0 + std::abs(kEps);
            const double dt = 0.4 * 0.6963 * h * h / (3.0 * D_max);

            std::vector<Snapshot> out;
            out.push_back({0, hd::computeMass(psi), hd::computeL2Squared(psi),
                           hd::computeMin(psi), hd::computeMax(psi)});

            int step = 0;
            for (int seg = 0; seg < kNumSegments; ++seg) {
                gridcalc::solve::diffuse(psi, D, dt, kStepsPerSegment, gridcalc::solve::RK4{});
                step += kStepsPerSegment;
                out.push_back({step, hd::computeMass(psi), hd::computeL2Squared(psi),
                               hd::computeMin(psi), hd::computeMax(psi)});
            }
            return out;
        }();
        return snaps;
    }
};

TEST_F(HeterogeneousDiffusionRun, MassConservation) {
    const auto& snaps = Snaps();
    ASSERT_EQ(snaps.size(), static_cast<std::size_t>(kNumSegments + 1));
    const double mass_initial = snaps.front().mass;
    for (std::size_t i = 1; i < snaps.size(); ++i) {
        EXPECT_NEAR(snaps[i].mass, mass_initial, 1e-10 * std::abs(mass_initial))
            << "mass drift at snapshot " << i;
    }
}

TEST_F(HeterogeneousDiffusionRun, L2NormMonotonicallyDecreases) {
    const auto& snaps = Snaps();
    // Strict inequality: each chunk integrates ~50 RK4 steps with non-zero
    // grad psi, so the L^2 entropy decays measurably.
    for (std::size_t i = 0; i + 1 < snaps.size(); ++i) {
        EXPECT_LT(snaps[i + 1].l2sq, snaps[i].l2sq)
            << "L^2 grew between snapshots " << i << " and " << (i + 1) << ": "
            << snaps[i].l2sq << " -> " << snaps[i + 1].l2sq;
    }
}

TEST_F(HeterogeneousDiffusionRun, PositivityPreserved) {
    const auto& snaps = Snaps();
    for (std::size_t i = 0; i < snaps.size(); ++i) {
        EXPECT_GE(snaps[i].psi_min, -1e-12)
            << "psi went negative at snapshot " << i << ": min=" << snaps[i].psi_min;
    }
}

}  // namespace
