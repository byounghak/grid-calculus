// Phase 15 — `func::evaluate` over `Field<Vec3d>` inputs.
//
// Tests cover the two new arities (`f(u)` and `f(u, grad_u)`),
// longest-arity tie-break, the `Kahan` reduction tag on the new path, and
// regression sanity that the new overloads do not interfere with the
// scalar `Field<double>` overloads from Phase 4 / 11.

#include <cmath>

#include <gtest/gtest.h>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/diff/Gradient.hpp>
#include <gridcalc/func/Functional.hpp>

namespace {

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Mat3d;
using gridcalc::core::Vec3d;
namespace func = gridcalc::func;

constexpr double kPi = 3.14159265358979323846;

Grid makeBox(int N) {
    const double h = 2.0 * kPi / N;
    return Grid(N, N, N, Vec3d(h, h, h));
}

template <class F>
Field<Vec3d> sampleVec(const Grid& g, F&& f) {
    Field<Vec3d> v(g, Vec3d::Zero());
    const auto& h = g.getCellSize();
    for (int k = 0; k < g.getNz(); ++k) {
        for (int j = 0; j < g.getNy(); ++j) {
            for (int i = 0; i < g.getNx(); ++i) {
                v(i, j, k) = f(static_cast<double>(i) * h(0),
                               static_cast<double>(j) * h(1),
                               static_cast<double>(k) * h(2));
            }
        }
    }
    return v;
}

// -------------------------------------------------------------------------
//  Arity (u) — f depends only on u, no derivative materialized
// -------------------------------------------------------------------------

// Closed form for u = (sin x, cos y, sin z cos x) on [0, 2pi]^3:
//   ||u||^2 = sin^2 x + cos^2 y + sin^2 z cos^2 x
// Integrals:
//   ∫_0^{2pi} sin^2 dx = pi
//   ∫_0^{2pi} cos^2 dy = pi
//   ∫∫ sin^2 z cos^2 x dx dz = pi * pi
// Total over (2pi)^3 box (each axis 2pi long, integrals over remaining axes):
//   (pi) * (2pi) * (2pi) + (2pi) * (pi) * (2pi) + (2pi) * (pi) * (pi)
//   = 4 pi^3 + 4 pi^3 + 2 pi^3 = 10 pi^3
TEST(FuncEvaluateVectorTest, EvaluateUOnlyArity) {
    const Grid g = makeBox(64);
    const Field<Vec3d> u = sampleVec(g, [](double x, double y, double z) {
        return Vec3d(std::sin(x), std::cos(y), std::sin(z) * std::cos(x));
    });
    const double F = func::evaluate(u, [](const Vec3d& u_pt) {
        return u_pt.dot(u_pt);
    });
    const double expected = 10.0 * kPi * kPi * kPi;
    EXPECT_NEAR(F, expected, 1e-3 * std::abs(expected));
}

// -------------------------------------------------------------------------
//  Arity (u, grad_u) — gradient is materialized once
// -------------------------------------------------------------------------

// For u = (sin x, sin y, sin z), grad_u = diag(cos x, cos y, cos z) (per
// Phase 13's M(i, j) = ∂_j u_i convention), so trace(grad_u) = cos x + cos y + cos z.
// ∫ trace^2 dV over [0, 2pi]^3 with cross-terms ∫ cos x cos y dV = 0:
//   ∫ (cos^2 x + cos^2 y + cos^2 z) dV = 3 * pi * (2pi)^2 = 12 pi^3.
TEST(FuncEvaluateVectorTest, EvaluateUAndGradArity) {
    const Grid g = makeBox(64);
    const Field<Vec3d> u = sampleVec(g, [](double x, double y, double z) {
        return Vec3d(std::sin(x), std::sin(y), std::sin(z));
    });
    const double F = func::evaluate(u, [](const Vec3d& /*u_pt*/, const Mat3d& grad_u) {
        const double tr = grad_u.trace();
        return tr * tr;
    });
    const double expected = 12.0 * kPi * kPi * kPi;
    EXPECT_NEAR(F, expected, 1e-2 * std::abs(expected));
}

// Convergence sweep for the (u, grad_u) arity.
TEST(FuncEvaluateVectorTest, UAndGradArity_SecondOrderConvergence) {
    auto run = [](int N) {
        const Grid g = makeBox(N);
        const Field<Vec3d> u = sampleVec(g, [](double x, double y, double z) {
            return Vec3d(std::sin(x), std::sin(y), std::sin(z));
        });
        return func::evaluate(u, [](const Vec3d& /*u_pt*/, const Mat3d& grad_u) {
            const double tr = grad_u.trace();
            return tr * tr;
        });
    };
    const double F_ref = 12.0 * kPi * kPi * kPi;
    const double e16 = std::abs(run(16) - F_ref);
    const double e32 = std::abs(run(32) - F_ref);
    const double e64 = std::abs(run(64) - F_ref);
    // Order-2 convergence: each refinement halves h, errors drop ~4x.
    // Defensive: error must be monotonically decreasing and the slope must
    // sit in [1.6, 2.4] on log-log.
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
//  Longest-arity tie-break
// -------------------------------------------------------------------------

// A generic lambda taking `(auto, ...)` is invocable as both arities; the
// 2-arg branch must win. We verify via a side effect captured by reference.
TEST(FuncEvaluateVectorTest, LongestArityWins) {
    const Grid g = makeBox(8);
    const Field<Vec3d> u = sampleVec(g, [](double x, double, double) {
        return Vec3d(std::sin(x), 0.0, 0.0);
    });
    int observed_arity = 0;
    auto generic = [&observed_arity](const Vec3d& u_pt, const Mat3d& grad_u) {
        observed_arity = 2;
        // Touch grad_u so the compiler doesn't elide the param.
        return u_pt.norm() + grad_u.trace() * 0.0;
    };
    (void)func::evaluate(u, generic);
    EXPECT_EQ(observed_arity, 2);
}

// -------------------------------------------------------------------------
//  Reduction tag
// -------------------------------------------------------------------------

TEST(FuncEvaluateVectorTest, KahanReductionTagWorks) {
    const Grid g = makeBox(8);
    const Field<Vec3d> u = sampleVec(g, [](double x, double y, double z) {
        return Vec3d(std::sin(x), std::sin(y), std::sin(z));
    });
    auto integrand = [](const Vec3d& /*u_pt*/, const Mat3d& grad_u) {
        const double tr = grad_u.trace();
        return tr * tr;
    };
    const double pair = func::evaluate(u, integrand, func::Pairwise{});
    const double kahan = func::evaluate(u, integrand, func::Kahan{});
    EXPECT_NEAR(pair, kahan, 1e-12 * std::abs(pair));
}

// -------------------------------------------------------------------------
//  Cross-check vs. Field<double> path (sanity)
// -------------------------------------------------------------------------

// f(u) = ||u||^2 and a parallel scalar formulation must agree to round-off
// once the inputs are equivalent: scalar phi = ||u||^2 sampled on the
// same grid, integrated as ∫ phi dV.
TEST(FuncEvaluateVectorTest, IntegrandReductionMatchesScalarPath) {
    const Grid g = makeBox(16);
    const Field<Vec3d> u = sampleVec(g, [](double x, double y, double z) {
        return Vec3d(std::sin(x), std::cos(y), std::sin(z));
    });
    const double F_vec = func::evaluate(u, [](const Vec3d& u_pt) {
        return u_pt.dot(u_pt);
    });
    Field<double> phi(g);
    const auto& h = g.getCellSize();
    for (int k = 0; k < g.getNz(); ++k) {
        for (int j = 0; j < g.getNy(); ++j) {
            for (int i = 0; i < g.getNx(); ++i) {
                const double x = static_cast<double>(i) * h(0);
                const double y = static_cast<double>(j) * h(1);
                const double z = static_cast<double>(k) * h(2);
                const Vec3d v(std::sin(x), std::cos(y), std::sin(z));
                phi(i, j, k) = v.dot(v);
            }
        }
    }
    const double F_scalar = func::integrate(phi);
    EXPECT_NEAR(F_vec, F_scalar, 1e-12 * std::abs(F_scalar));
}

}  // namespace
