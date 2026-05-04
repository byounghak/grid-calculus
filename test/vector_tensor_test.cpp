// Phase 13 — vector and tensor field operators.
//
// Layout:
//   * `MatGradient*`  — rank-2 gradient `Field<Vec3d> -> Field<Mat3d>`.
//   * `MatDivergence*` — rank-2 divergence `Field<Mat3d> -> Field<Vec3d>`.
//   * `Trace*`, `SingleContract*`, `DoubleContract*` — contraction primitives.
//   * `VectorIdentity*` — the three roadmap-acceptance identities.
//
// Group 3 covers the `MatGradient*` block; later groups append.

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/diff/Gradient.hpp>

namespace {

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Mat3d;
using gridcalc::core::Vec3d;
using gridcalc::diff::gradient;

constexpr double kPi = 3.14159265358979323846;

// Builds a `Field<Vec3d>` by sampling a vector-valued callable on the
// grid in i-fastest order. (Matches the Phase 1 `Field(grid, f)`
// constructor's loop ordering.)
template <class F>
Field<Vec3d> sampleVectorField(const Grid& g, F&& f) {
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

// =============================================================================
// MatGradient
// =============================================================================

TEST(MatGradientTest, OnLinearField) {
    // v = (a*x, b*y, c*z) -> M = diag(a, b, c) at every interior point.
    // Linear input + 2nd-order central differences -> exact (round-off only).
    constexpr int N = 8;
    const double L = 2.0 * kPi;
    const double h = L / N;
    Grid g(N, N, N, Vec3d(h, h, h));

    const double a = 1.5, b = -0.5, c = 2.25;
    auto v = sampleVectorField(g, [&](double x, double y, double z) {
        return Vec3d(a * x, b * y, c * z);
    });
    auto M = gradient(v);

    // Linear-field derivatives are exact under the central stencil only at
    // interior points (boundary wrapping introduces a finite-step bias on a
    // linear-but-non-periodic input). Check a few interior cells.
    for (int k = 1; k + 1 < N; ++k) {
        for (int j = 1; j + 1 < N; ++j) {
            for (int i = 1; i + 1 < N; ++i) {
                const Mat3d& Mij = M(i, j, k);
                EXPECT_NEAR(Mij(0, 0), a, 1e-12) << "i,j,k=" << i << "," << j << "," << k;
                EXPECT_NEAR(Mij(1, 1), b, 1e-12);
                EXPECT_NEAR(Mij(2, 2), c, 1e-12);
                EXPECT_NEAR(Mij(0, 1), 0.0, 1e-12);
                EXPECT_NEAR(Mij(0, 2), 0.0, 1e-12);
                EXPECT_NEAR(Mij(1, 0), 0.0, 1e-12);
                EXPECT_NEAR(Mij(1, 2), 0.0, 1e-12);
                EXPECT_NEAR(Mij(2, 0), 0.0, 1e-12);
                EXPECT_NEAR(Mij(2, 1), 0.0, 1e-12);
            }
        }
    }
}

TEST(MatGradientTest, IndexConvention) {
    // v = (sin(y), 0, 0). Then ∂_y v_x = cos(y); all other components are 0.
    // Under the M(i,j) = ∂_j v_i convention, the only non-zero entry is
    // M(0, 1) = cos(y).
    constexpr int N = 32;
    const double L = 2.0 * kPi;
    const double h = L / N;
    Grid g(N, N, N, Vec3d(h, h, h));
    auto v = sampleVectorField(g, [](double, double y, double) {
        return Vec3d(std::sin(y), 0.0, 0.0);
    });
    auto M = gradient(v);

    // Sample one interior point; the FD error is O(h^2) in cos(y), well
    // below the 1e-2 sanity tolerance below.
    const int i0 = 7, j0 = 11, k0 = 5;
    const double y = static_cast<double>(j0) * h;
    const Mat3d& Mij = M(i0, j0, k0);
    EXPECT_NEAR(Mij(0, 1), std::cos(y), 1e-2);
    // Every other entry should be at machine precision (the Vec3d input
    // has identically-zero v_y / v_z components).
    EXPECT_NEAR(Mij(0, 0), 0.0, 1e-12);
    EXPECT_NEAR(Mij(0, 2), 0.0, 1e-12);
    EXPECT_NEAR(Mij(1, 0), 0.0, 1e-12);
    EXPECT_NEAR(Mij(1, 1), 0.0, 1e-12);
    EXPECT_NEAR(Mij(1, 2), 0.0, 1e-12);
    EXPECT_NEAR(Mij(2, 0), 0.0, 1e-12);
    EXPECT_NEAR(Mij(2, 1), 0.0, 1e-12);
    EXPECT_NEAR(Mij(2, 2), 0.0, 1e-12);
}

// Helper: max-norm of (computed gradient - analytical gradient) for the
// trig manufactured solution
//   v(x, y, z) = (sin(x), sin(y), sin(z))
// whose analytical rank-2 gradient is M_ref = diag(cos x, cos y, cos z).
double matGradientErrorAtN(int N) {
    const double L = 2.0 * kPi;
    const double h = L / N;
    Grid g(N, N, N, Vec3d(h, h, h));
    auto v = sampleVectorField(g, [](double x, double y, double z) {
        return Vec3d(std::sin(x), std::sin(y), std::sin(z));
    });
    auto M = gradient(v);

    double err = 0.0;
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                const double x = static_cast<double>(i) * h;
                const double y = static_cast<double>(j) * h;
                const double z = static_cast<double>(k) * h;
                Mat3d M_ref = Mat3d::Zero();
                M_ref(0, 0) = std::cos(x);
                M_ref(1, 1) = std::cos(y);
                M_ref(2, 2) = std::cos(z);
                err = std::max(err, (M(i, j, k) - M_ref).cwiseAbs().maxCoeff());
            }
        }
    }
    return err;
}

TEST(MatGradientTest, ConvergenceOrder2) {
    const double e16 = matGradientErrorAtN(16);
    const double e32 = matGradientErrorAtN(32);
    const double e64 = matGradientErrorAtN(64);

    const double slope_16_32 = std::log2(e16 / e32);
    const double slope_32_64 = std::log2(e32 / e64);
    EXPECT_GT(slope_16_32, 1.7);
    EXPECT_LT(slope_16_32, 2.3);
    EXPECT_GT(slope_32_64, 1.7);
    EXPECT_LT(slope_32_64, 2.3);
}

TEST(MatGradientTest, Order4IsTighterThanOrder2) {
    // Same trig manufactured solution; switch to the Phase 7 4th-order
    // stencil and confirm the error at fixed N drops by at least a factor
    // of 5. Loose threshold so we are not racing the inversely-scaled
    // round-off floor.
    constexpr int N = 32;
    const double L = 2.0 * kPi;
    const double h = L / N;
    Grid g(N, N, N, Vec3d(h, h, h));
    auto v = sampleVectorField(g, [](double x, double y, double z) {
        return Vec3d(std::sin(x), std::sin(y), std::sin(z));
    });
    auto M2 = gradient<2>(v);
    auto M4 = gradient<4>(v);

    double err2 = 0.0, err4 = 0.0;
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                const double x = static_cast<double>(i) * h;
                const double y = static_cast<double>(j) * h;
                const double z = static_cast<double>(k) * h;
                Mat3d M_ref = Mat3d::Zero();
                M_ref(0, 0) = std::cos(x);
                M_ref(1, 1) = std::cos(y);
                M_ref(2, 2) = std::cos(z);
                err2 = std::max(err2, (M2(i, j, k) - M_ref).cwiseAbs().maxCoeff());
                err4 = std::max(err4, (M4(i, j, k) - M_ref).cwiseAbs().maxCoeff());
            }
        }
    }
    EXPECT_LT(err4, err2 / 5.0);
}

}  // namespace
