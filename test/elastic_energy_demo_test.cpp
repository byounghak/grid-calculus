// Phase 15 demo CI gate.
//
// Mirrors the Phase 12 / Phase 14 pattern: runs a stripped-down
// version of `examples/elastic_energy.cpp` (just the helpers from
// `examples/common/elastic_energy.hpp` + `func::evaluate`, no VTK
// output, no CLI parsing) at a small grid extent on every PR. Asserts
// the relative error is within a loose bound at the small grid size.

#include <cmath>

#include <gtest/gtest.h>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/func/Functional.hpp>

#include <common/elastic_energy.hpp>

namespace {

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Mat3d;
using gridcalc::core::Vec3d;
namespace ee = gridcalc::examples::elastic_energy;
namespace func = gridcalc::func;

constexpr double kPi = 3.14159265358979323846;

TEST(ElasticEnergyDemoTest, ComputesEnergyToWithinPercentAtN24) {
    constexpr int N = 24;
    const double h = 2.0 * kPi / N;
    Grid grid(N, N, N, Vec3d(h, h, h));

    const auto lame = ee::lameFromYoungPoisson(/*E=*/1.0, /*nu=*/0.3);
    const double lambda = lame.first;
    const double mu = lame.second;
    const double alpha = 0.01;

    const Field<Vec3d> u = ee::makeUniaxialPeriodicDisplacement(grid, alpha, /*axis=*/0);
    const double F_h = func::evaluate(u, [lambda, mu](const Vec3d&, const Mat3d& gu) {
        return ee::linearElasticEnergyDensity(lambda, mu, gu);
    });
    const double F_ref =
        ee::analyticalLinearElasticEnergy(grid, alpha, lambda, mu, /*axis=*/0);
    const double rel = std::abs(F_h - F_ref) / std::abs(F_ref);
    // At N = 24 the order-2 leading constant gives ~2.3 x 10^-2 relative
    // error; gate at 5 x 10^-2 to leave headroom for floating-point noise
    // across compilers.
    EXPECT_LT(rel, 5e-2) << "F_h = " << F_h << ", F_ref = " << F_ref << ", rel = " << rel;
    EXPECT_GT(F_h, 0.0);  // sanity: energy is non-trivial and positive.
}

}  // namespace
