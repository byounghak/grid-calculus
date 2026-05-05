// Phase 15 — linear-elastic-energy acceptance test.
//
// Discharges the roadmap acceptance bar: the elastic energy of a
// (periodic) "stretched cubic block" matches the analytical value. The
// rigid-stretch geometry is impossible on the only-`IndexPolicy::Periodic`
// surface; the periodic substitute is a uniaxial dilation wave
// `u_axis = alpha sin(k0 * x_axis)` whose closed-form energy is
// `F_ref = 2 pi^3 (lambda + 2 mu) alpha^2` on the standard
// `[0, 2 pi]^3` box (see `examples/common/elastic_energy.hpp`).

#include <cmath>

#include <gtest/gtest.h>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/diff/Gradient.hpp>
#include <gridcalc/func/Functional.hpp>
#include <gridcalc/func/Integrate.hpp>
#include <gridcalc/tensor/Expressions.hpp>

#include <common/elastic_energy.hpp>

namespace {

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Mat3d;
using gridcalc::core::Vec3d;
namespace func = gridcalc::func;
namespace expr = gridcalc::tensor::expr;
namespace ee = gridcalc::examples::elastic_energy;

constexpr double kPi = 3.14159265358979323846;
constexpr double kE = 1.0;
constexpr double kNu = 0.3;
constexpr double kAlpha = 0.01;

Grid makeStandardBox(int N) {
    const double h = 2.0 * kPi / N;
    return Grid(N, N, N, Vec3d(h, h, h));
}

double computeElasticEnergy(const Grid& g, double alpha, double lambda, double mu, int axis) {
    const Field<Vec3d> u = ee::makeUniaxialPeriodicDisplacement(g, alpha, axis);
    return func::evaluate(u, [lambda, mu](const Vec3d&, const Mat3d& grad_u) {
        return ee::linearElasticEnergyDensity(lambda, mu, grad_u);
    });
}

// -------------------------------------------------------------------------
//  Pointwise sanity (no grid)
// -------------------------------------------------------------------------

TEST(ElasticEnergyTest, EnergyDensityFormulaAtZeroStrain) {
    const auto lame = ee::lameFromYoungPoisson(kE, kNu);
    const double lambda = lame.first;
    const double mu = lame.second;
    EXPECT_DOUBLE_EQ(ee::linearElasticEnergyDensity(lambda, mu, Mat3d::Zero()), 0.0);
}

TEST(ElasticEnergyTest, EnergyDensityFormulaAtUnitStrainAlongX) {
    // grad_u = diag(1, 0, 0) -> eps = diag(1, 0, 0); trace = 1, eps:eps = 1.
    // energy density = 0.5 * (lambda * 1 + 2 mu * 1) = 0.5 * (lambda + 2 mu).
    const auto lame = ee::lameFromYoungPoisson(kE, kNu);
    const double lambda = lame.first;
    const double mu = lame.second;
    Mat3d gu = Mat3d::Zero();
    gu(0, 0) = 1.0;
    const double expected = 0.5 * (lambda + 2.0 * mu);
    EXPECT_NEAR(ee::linearElasticEnergyDensity(lambda, mu, gu), expected, 1e-15);
}

TEST(ElasticEnergyTest, LameRoundTrip) {
    const auto lame = ee::lameFromYoungPoisson(kE, kNu);
    const double lambda = lame.first;
    const double mu = lame.second;
    // For E = 1, nu = 0.3: lambda = 0.3 / (1.3 * 0.4) = 0.5769230769...
    //                     mu     = 1   / (2 * 1.3)  = 0.3846153846...
    EXPECT_NEAR(lambda, 0.3 / (1.3 * 0.4), 1e-12);
    EXPECT_NEAR(mu, 1.0 / (2.0 * 1.3), 1e-12);
}

TEST(ElasticEnergyTest, LameRejectsBadInputs) {
    EXPECT_THROW((void)ee::lameFromYoungPoisson(0.0, 0.3), std::invalid_argument);
    EXPECT_THROW((void)ee::lameFromYoungPoisson(-1.0, 0.3), std::invalid_argument);
    EXPECT_THROW((void)ee::lameFromYoungPoisson(1.0, 0.5), std::invalid_argument);
    EXPECT_THROW((void)ee::lameFromYoungPoisson(1.0, -1.0), std::invalid_argument);
}

// -------------------------------------------------------------------------
//  Acceptance — closed-form match at N = 64
// -------------------------------------------------------------------------

TEST(ElasticEnergyTest, MatchesAnalyticalAtN64) {
    // Order-2 central differences leave a leading constant ~0.3 in the
    // relative-error series; at h = 2 pi / 64 the discrete energy lands
    // ~3 x 10^-3 from the closed form. The convergence sweep below
    // verifies the rate; this test verifies the absolute value at the
    // recommended demo grid size with a 5 x 10^-3 safety margin.
    const Grid g = makeStandardBox(64);
    const auto lame = ee::lameFromYoungPoisson(kE, kNu);
    const double lambda = lame.first;
    const double mu = lame.second;
    const double F_h = computeElasticEnergy(g, kAlpha, lambda, mu, /*axis=*/0);
    const double F_ref = ee::analyticalLinearElasticEnergy(g, kAlpha, lambda, mu, /*axis=*/0);
    const double rel = std::abs(F_h - F_ref) / std::abs(F_ref);
    EXPECT_LT(rel, 5e-3) << "F_h = " << F_h << ", F_ref = " << F_ref << ", rel = " << rel;
}

// -------------------------------------------------------------------------
//  Convergence sweep — order 2
// -------------------------------------------------------------------------

TEST(ElasticEnergyTest, SecondOrderConvergence) {
    const auto lame = ee::lameFromYoungPoisson(kE, kNu);
    const double lambda = lame.first;
    const double mu = lame.second;
    auto err = [&](int N) {
        const Grid g = makeStandardBox(N);
        const double F_h = computeElasticEnergy(g, kAlpha, lambda, mu, /*axis=*/0);
        const double F_ref = ee::analyticalLinearElasticEnergy(g, kAlpha, lambda, mu, /*axis=*/0);
        return std::abs(F_h - F_ref);
    };
    const double e16 = err(16);
    const double e32 = err(32);
    const double e64 = err(64);
    EXPECT_GT(e16, e32);
    EXPECT_GT(e32, e64);
    const double slope_low = std::log2(e16 / e32);
    const double slope_high = std::log2(e32 / e64);
    EXPECT_GE(slope_low, 1.6);
    EXPECT_LE(slope_low, 2.4);
    EXPECT_GE(slope_high, 1.6);
    EXPECT_LE(slope_high, 2.4);
}

// -------------------------------------------------------------------------
//  Symmetry — energy is axis-independent
// -------------------------------------------------------------------------

TEST(ElasticEnergyTest, EachAxisGivesSameAnalyticalEnergy) {
    const Grid g = makeStandardBox(32);
    const auto lame = ee::lameFromYoungPoisson(kE, kNu);
    const double lambda = lame.first;
    const double mu = lame.second;
    const double Fx = computeElasticEnergy(g, kAlpha, lambda, mu, /*axis=*/0);
    const double Fy = computeElasticEnergy(g, kAlpha, lambda, mu, /*axis=*/1);
    const double Fz = computeElasticEnergy(g, kAlpha, lambda, mu, /*axis=*/2);
    EXPECT_NEAR(Fy, Fx, 1e-12 * std::abs(Fx));
    EXPECT_NEAR(Fz, Fx, 1e-12 * std::abs(Fx));
}

// -------------------------------------------------------------------------
//  Cross-check via the contraction expression-template layer
// -------------------------------------------------------------------------

TEST(ElasticEnergyTest, MatchesViaExpressionLayer) {
    const Grid g = makeStandardBox(32);
    const auto lame = ee::lameFromYoungPoisson(kE, kNu);
    const double lambda = lame.first;
    const double mu = lame.second;

    const Field<Vec3d> u = ee::makeUniaxialPeriodicDisplacement(g, kAlpha, /*axis=*/0);
    const Field<Mat3d> grad_u = gridcalc::diff::gradient(u);

    // Path A — `func::evaluate` with the per-point Lamé closure.
    const double F_eval = func::evaluate(u, [lambda, mu](const Vec3d&, const Mat3d& gu) {
        return ee::linearElasticEnergyDensity(lambda, mu, gu);
    });

    // Path B — fused integrate over the contraction-ET expression.
    auto eps_expr = expr::sym(expr::field(grad_u));
    auto trace_eps = expr::trace(eps_expr);
    auto integrand = 0.5 * lambda * (trace_eps * trace_eps) +
                     mu * expr::doubleContract(eps_expr, eps_expr);
    const double F_expr = func::integrate(integrand);

    EXPECT_NEAR(F_eval, F_expr, 1e-12 * std::abs(F_eval));
}

}  // namespace
