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

}  // namespace
