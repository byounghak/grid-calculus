// Phase 14 — fvm::cellLaplacian tests.
//
// Layout (filled in across Group 3 / Group 4):
//   * Group 3: constant-D agreement gate vs. diff::laplacian + zero-RHS sanity.
//   * Group 4: harmonic-mean face-D + heterogeneous-D convergence + mass conservation.

#include <cmath>

#include <gtest/gtest.h>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/diff/Laplacian.hpp>
#include <gridcalc/fvm/CellLaplacian.hpp>

namespace {

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;

constexpr double kPi = 3.14159265358979323846;

Field<double> sampleScalarField(const Grid& g,
                                double (*f)(double, double, double)) {
    return Field<double>(g, f);
}

TEST(FvmCellLaplacianTest, AgreesWithFdLaplacianAtConstantUnitD) {
    constexpr int N = 16;
    const double L = 2.0 * kPi;
    const double h = L / N;
    Grid g(N, N, N, Vec3d(h, h, h));
    auto psi = sampleScalarField(g, [](double x, double y, double z) {
        return std::sin(x) * std::cos(y) * std::sin(z);
    });
    Field<double> D(g, 1.0);

    auto fvm_lap = gridcalc::fvm::cellLaplacian(psi, D);
    auto fd_lap = gridcalc::diff::laplacian(psi);

    double max_diff = 0.0;
    for (std::size_t idx = 0; idx < psi.getSize(); ++idx) {
        max_diff = std::max(max_diff, std::abs(fvm_lap.data()[idx] - fd_lap.data()[idx]));
    }
    EXPECT_LT(max_diff, 1e-12);
}

TEST(FvmCellLaplacianTest, AgreesWithFdLaplacianAtConstantD) {
    constexpr int N = 16;
    const double L = 2.0 * kPi;
    const double h = L / N;
    Grid g(N, N, N, Vec3d(h, h, h));
    auto psi = sampleScalarField(g, [](double x, double y, double z) {
        return std::sin(x) * std::cos(y) * std::sin(z);
    });
    const double c = 0.7;
    Field<double> D(g, c);

    auto fvm_lap = gridcalc::fvm::cellLaplacian(psi, D);
    auto fd_lap = gridcalc::diff::laplacian(psi);

    double max_diff = 0.0;
    for (std::size_t idx = 0; idx < psi.getSize(); ++idx) {
        max_diff = std::max(max_diff,
                            std::abs(fvm_lap.data()[idx] - c * fd_lap.data()[idx]));
    }
    EXPECT_LT(max_diff, 1e-12);
}

TEST(FvmCellLaplacianTest, ConstantPsiGivesZeroRhs) {
    Grid g(8, 8, 8, Vec3d(1.0, 1.0, 1.0));
    Field<double> psi(g, 5.0);
    // Non-trivial D; the RHS should still be zero everywhere because
    // the gradient of a constant field is zero, so every face flux is
    // zero.
    auto D = sampleScalarField(g, [](double x, double y, double z) {
        return 1.0 + 0.5 * std::sin(x) * std::cos(y) * std::sin(z);
    });

    auto rhs = gridcalc::fvm::cellLaplacian(psi, D);
    double max_abs = 0.0;
    for (std::size_t idx = 0; idx < rhs.getSize(); ++idx) {
        max_abs = std::max(max_abs, std::abs(rhs.data()[idx]));
    }
    EXPECT_LT(max_abs, 1e-12);
}

TEST(FvmCellLaplacianTest, HarmonicMeanFaceCoefficient) {
    // 4-cell strip along x with a D-jump and a known psi pattern.
    // D = [1, 4, 4, 1], psi = [0, 1, 1, 0], hx = 1, periodic.
    // For cell i = 1:
    //   D_xp = hmean(D(1), D(2)) = hmean(4, 4) = 4
    //   D_xm = hmean(D(0), D(1)) = hmean(1, 4) = 8/5 = 1.6
    //   rhs_x = (4 * (psi(2) - psi(1)) - 1.6 * (psi(1) - psi(0))) / hx^2
    //         = (4 * 0 - 1.6 * 1) / 1 = -1.6
    // For cell i = 0 (psi = 0, D_xm wraps to D(3) = 1, psi_xm wraps to psi(3) = 0):
    //   D_xp = hmean(D(0), D(1)) = hmean(1, 4) = 1.6
    //   D_xm = hmean(D(3), D(0)) = hmean(1, 1) = 1
    //   rhs_x = (1.6 * (psi(1) - psi(0)) - 1 * (psi(0) - psi(3))) / 1
    //         = (1.6 * 1 - 1 * 0) = 1.6
    Grid g(4, 1, 1, Vec3d(1.0, 1.0, 1.0));
    Field<double> psi(g, 0.0);
    Field<double> D(g, 0.0);
    psi(0, 0, 0) = 0.0;
    psi(1, 0, 0) = 1.0;
    psi(2, 0, 0) = 1.0;
    psi(3, 0, 0) = 0.0;
    D(0, 0, 0) = 1.0;
    D(1, 0, 0) = 4.0;
    D(2, 0, 0) = 4.0;
    D(3, 0, 0) = 1.0;

    auto rhs = gridcalc::fvm::cellLaplacian(psi, D);
    EXPECT_DOUBLE_EQ(rhs(0, 0, 0), 1.6);
    EXPECT_DOUBLE_EQ(rhs(1, 0, 0), -1.6);
    EXPECT_DOUBLE_EQ(rhs(2, 0, 0), -1.6);
    EXPECT_DOUBLE_EQ(rhs(3, 0, 0), 1.6);
}

// Manufactured solution for the 1D-flavored heterogeneous-D operator:
//   D(x, y, z) = D0 + eps * cos(x)
//   psi(x, y, z) = sin(x)
// The cell-Laplacian collapses to the x-axis term:
//   div(D grad psi) = d_x [D * d_x psi] = d_x [(D0 + eps cos x) * cos x]
//                   = d_x [D0 cos x + eps cos^2 x]
//                   = -D0 sin x - 2 eps cos x sin x
//                   = -D0 sin x - eps sin(2x).
// y and z derivatives are zero because psi has no y, z dependence and
// D's y, z dependence is also zero in the manufactured solution.
constexpr double kD0 = 1.0;
constexpr double kEps = 0.5;

double manufacturedRhs(double x) {
    return -kD0 * std::sin(x) - kEps * std::sin(2.0 * x);
}

double heterogeneousErrorAtN(int N) {
    const double L = 2.0 * kPi;
    const double h = L / N;
    Grid g(N, N, N, Vec3d(h, h, h));
    auto psi = sampleScalarField(g, [](double x, double, double) {
        return std::sin(x);
    });
    auto D = sampleScalarField(g, [](double x, double, double) {
        return kD0 + kEps * std::cos(x);
    });
    auto rhs = gridcalc::fvm::cellLaplacian(psi, D);

    double err = 0.0;
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                const double x = static_cast<double>(i) * h;
                err = std::max(err, std::abs(rhs(i, j, k) - manufacturedRhs(x)));
            }
        }
    }
    return err;
}

TEST(FvmCellLaplacianTest, AnalyticalSolutionConvergence_Order2) {
    const double e16 = heterogeneousErrorAtN(16);
    const double e32 = heterogeneousErrorAtN(32);
    const double e64 = heterogeneousErrorAtN(64);

    const double slope_16_32 = std::log2(e16 / e32);
    const double slope_32_64 = std::log2(e32 / e64);
    EXPECT_GT(slope_16_32, 1.7);
    EXPECT_LT(slope_16_32, 2.3);
    EXPECT_GT(slope_32_64, 1.7);
    EXPECT_LT(slope_32_64, 2.3);
}

TEST(FvmCellLaplacianTest, MassConservation_PeriodicSum) {
    // Sum_{i,j,k} fvm::cellLaplacian(psi, D)(i, j, k) = 0 to round-off
    // on a periodic domain (telescoping fluxes — every face flux is
    // counted once positive and once negative across the two adjacent
    // cells). Holds for arbitrary psi and arbitrary D > 0.
    constexpr int N = 32;
    const double L = 2.0 * kPi;
    const double h = L / N;
    Grid g(N, N, N, Vec3d(h, h, h));
    auto psi = sampleScalarField(g, [](double x, double y, double z) {
        return std::sin(2.0 * x) * std::cos(y) + 0.3 * std::sin(z);
    });
    auto D = sampleScalarField(g, [](double x, double y, double z) {
        return 1.0 + 0.7 * (std::sin(x) * std::cos(y) * std::sin(z));
    });

    auto rhs = gridcalc::fvm::cellLaplacian(psi, D);
    double sum = 0.0;
    double max_abs = 0.0;
    for (std::size_t idx = 0; idx < rhs.getSize(); ++idx) {
        sum += rhs.data()[idx];
        max_abs = std::max(max_abs, std::abs(rhs.data()[idx]));
    }
    // Tolerance scaled by point count and per-point magnitude (flat
    // floating-point sum accumulates O(N^3) round-off).
    EXPECT_LT(std::abs(sum), 1e-10 * max_abs * static_cast<double>(rhs.getSize()));
}

}  // namespace
