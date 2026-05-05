// Phase 15 review-fix #1: Grid-equality precondition on the eager
// contraction primitives `tensor::singleContract` and
// `tensor::doubleContract`. Same defensive-contract gap that PR #20
// closed for `solve::integrate` / `axpyFresh` and PR #23 closed for
// `fvm::cellLaplacian` / `solve::diffuse`, just left open in Phase 13
// until the Phase-15 review surfaced it.
//
// Before 0.15.0 these primitives looped `idx < a.getGrid().getSize()`
// reading `b.data()[idx]` with no check; equal-`getSize()` /
// different-extent grids paired physically unrelated cells, and a
// smaller `b` triggered an OOB read.

#include <gtest/gtest.h>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/tensor/Contraction.hpp>
#include <gridcalc/tensor/Expressions.hpp>

namespace {

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Mat3d;
using gridcalc::core::Vec3d;
namespace tensor = gridcalc::tensor;
namespace expr = gridcalc::tensor::expr;

constexpr double kPi = 3.14159265358979323846;

Grid makeBox(int Nx, int Ny, int Nz, double hx, double hy, double hz) {
    return Grid(Nx, Ny, Nz, Vec3d(hx, hy, hz));
}

Grid makeStandardBox(int N) {
    const double h = 2.0 * kPi / N;
    return makeBox(N, N, N, h, h, h);
}

// --- Three flavors of mismatch -------------------------------------------
// 1. Smaller storage on the rhs: `singleContract`'s flat loop would OOB.
// 2. Same total size, different extents: would silently pair unrelated cells.
// 3. Same shape, different cell size: would silently use wrong physical h.

TEST(TensorGridMismatchTest, SingleContract_RejectsSmallerRhs) {
    Field<Mat3d> a(makeStandardBox(8), Mat3d::Identity());
    Field<Mat3d> b(makeStandardBox(4), Mat3d::Identity());
    EXPECT_THROW((void)tensor::singleContract(a, b), std::invalid_argument);
}

TEST(TensorGridMismatchTest, SingleContract_RejectsSameSizeDifferentExtents) {
    // 8*8*8 == 4*16*8 == 512 cells, but the per-axis extents differ.
    const Grid g_left = makeBox(8, 8, 8, 0.1, 0.1, 0.1);
    const Grid g_right = makeBox(4, 16, 8, 0.1, 0.1, 0.1);
    Field<Mat3d> a(g_left, Mat3d::Identity());
    Field<Mat3d> b(g_right, Mat3d::Identity());
    EXPECT_THROW((void)tensor::singleContract(a, b), std::invalid_argument);
}

TEST(TensorGridMismatchTest, SingleContract_RejectsSameShapeDifferentSpacing) {
    Field<Mat3d> a(makeBox(8, 8, 8, 0.1, 0.1, 0.1), Mat3d::Identity());
    Field<Mat3d> b(makeBox(8, 8, 8, 0.2, 0.1, 0.1), Mat3d::Identity());
    EXPECT_THROW((void)tensor::singleContract(a, b), std::invalid_argument);
}

TEST(TensorGridMismatchTest, DoubleContract_RejectsSmallerRhs) {
    Field<Mat3d> a(makeStandardBox(8), Mat3d::Identity());
    Field<Mat3d> b(makeStandardBox(4), Mat3d::Identity());
    EXPECT_THROW((void)tensor::doubleContract(a, b), std::invalid_argument);
}

TEST(TensorGridMismatchTest, DoubleContract_RejectsSameSizeDifferentExtents) {
    const Grid g_left = makeBox(8, 8, 8, 0.1, 0.1, 0.1);
    const Grid g_right = makeBox(4, 16, 8, 0.1, 0.1, 0.1);
    Field<Mat3d> a(g_left, Mat3d::Identity());
    Field<Mat3d> b(g_right, Mat3d::Identity());
    EXPECT_THROW((void)tensor::doubleContract(a, b), std::invalid_argument);
}

TEST(TensorGridMismatchTest, DoubleContract_RejectsSameShapeDifferentSpacing) {
    Field<Mat3d> a(makeBox(8, 8, 8, 0.1, 0.1, 0.1), Mat3d::Identity());
    Field<Mat3d> b(makeBox(8, 8, 8, 0.1, 0.2, 0.1), Mat3d::Identity());
    EXPECT_THROW((void)tensor::doubleContract(a, b), std::invalid_argument);
}

// --- Happy path still works ---------------------------------------------

TEST(TensorGridMismatchTest, MatchingGridsStillWork) {
    Field<Mat3d> a(makeStandardBox(4), Mat3d::Identity());
    Field<Mat3d> b(makeStandardBox(4), Mat3d::Identity());
    Field<Mat3d> single = tensor::singleContract(a, b);
    Field<double> double_c = tensor::doubleContract(a, b);
    // Identity * Identity = Identity (per cell); Identity : Identity = 3.
    for (int k = 0; k < a.getGrid().getNz(); ++k) {
        for (int j = 0; j < a.getGrid().getNy(); ++j) {
            for (int i = 0; i < a.getGrid().getNx(); ++i) {
                EXPECT_TRUE(single(i, j, k).isApprox(Mat3d::Identity(), 1e-15));
                EXPECT_DOUBLE_EQ(double_c(i, j, k), 3.0);
            }
        }
    }
}

// --- Error message includes context tag ----------------------------------

// --- ET binary nodes (review fix #2) -------------------------------------
// Without grid validation, Plus / Minus / Mul / SingleContract /
// DoubleContract would silently mis-pair cells via the periodic-wrap
// policy on `Field::operator()` and integrate the wrong physics.

TEST(TensorGridMismatchTest, ETPlus_RejectsSmallerRhs) {
    Field<Mat3d> a(makeStandardBox(8), Mat3d::Identity());
    Field<Mat3d> b(makeStandardBox(4), Mat3d::Identity());
    auto la = expr::field(a);
    auto lb = expr::field(b);
    EXPECT_THROW((void)(la + lb), std::invalid_argument);
}

TEST(TensorGridMismatchTest, ETMinus_RejectsSameSizeDifferentExtents) {
    const Grid g_left = makeBox(8, 8, 8, 0.1, 0.1, 0.1);
    const Grid g_right = makeBox(4, 16, 8, 0.1, 0.1, 0.1);
    Field<Mat3d> a(g_left, Mat3d::Identity());
    Field<Mat3d> b(g_right, Mat3d::Identity());
    auto la = expr::field(a);
    auto lb = expr::field(b);
    EXPECT_THROW((void)(la - lb), std::invalid_argument);
}

TEST(TensorGridMismatchTest, ETMul_RejectsSameShapeDifferentSpacing) {
    Field<Mat3d> a(makeBox(8, 8, 8, 0.1, 0.1, 0.1), Mat3d::Identity());
    Field<Mat3d> b(makeBox(8, 8, 8, 0.2, 0.1, 0.1), Mat3d::Identity());
    auto trA = expr::trace(expr::field(a));
    auto trB = expr::trace(expr::field(b));
    EXPECT_THROW((void)(trA * trB), std::invalid_argument);
}

TEST(TensorGridMismatchTest, ETSingleContract_RejectsSmallerRhs) {
    Field<Mat3d> a(makeStandardBox(8), Mat3d::Identity());
    Field<Mat3d> b(makeStandardBox(4), Mat3d::Identity());
    auto la = expr::field(a);
    auto lb = expr::field(b);
    EXPECT_THROW((void)expr::singleContract(la, lb), std::invalid_argument);
}

TEST(TensorGridMismatchTest, ETDoubleContract_RejectsSameSizeDifferentExtents) {
    const Grid g_left = makeBox(8, 8, 8, 0.1, 0.1, 0.1);
    const Grid g_right = makeBox(4, 16, 8, 0.1, 0.1, 0.1);
    Field<Mat3d> a(g_left, Mat3d::Identity());
    Field<Mat3d> b(g_right, Mat3d::Identity());
    auto la = expr::field(a);
    auto lb = expr::field(b);
    EXPECT_THROW((void)expr::doubleContract(la, lb), std::invalid_argument);
}

TEST(TensorGridMismatchTest, ETDoubleContract_RejectsSameShapeDifferentSpacing) {
    Field<Mat3d> a(makeBox(8, 8, 8, 0.1, 0.1, 0.1), Mat3d::Identity());
    Field<Mat3d> b(makeBox(8, 8, 8, 0.1, 0.2, 0.1), Mat3d::Identity());
    auto la = expr::field(a);
    auto lb = expr::field(b);
    EXPECT_THROW((void)expr::doubleContract(la, lb), std::invalid_argument);
}

TEST(TensorGridMismatchTest, ETPlusWithIdentityField_RejectsMismatch) {
    // IdentityField + Leaf<Mat3d>: both expose grid(); the binary node
    // must validate them.
    Field<Mat3d> a(makeStandardBox(8), Mat3d::Identity());
    const Grid g_other = makeStandardBox(4);
    auto leaf = expr::field(a);
    auto idf = expr::identityField(g_other);
    EXPECT_THROW((void)(leaf + idf), std::invalid_argument);
}

TEST(TensorGridMismatchTest, ETBinary_HappyPathStillWorks) {
    Field<Mat3d> a(makeStandardBox(4), Mat3d::Identity());
    Field<Mat3d> b(makeStandardBox(4), Mat3d::Identity());
    auto la = expr::field(a);
    auto lb = expr::field(b);
    // Build all five binary nodes; they must all construct without throw.
    EXPECT_NO_THROW({
        auto sum = la + lb;
        auto diff = la - lb;
        auto sing = expr::singleContract(la, lb);
        auto dub = expr::doubleContract(la, lb);
        auto trA = expr::trace(la);
        auto trB = expr::trace(lb);
        auto prod = trA * trB;
        (void)sum; (void)diff; (void)sing; (void)dub; (void)prod;
    });
}

TEST(TensorGridMismatchTest, ETErrorMessageNamesNode) {
    Field<Mat3d> a(makeStandardBox(8), Mat3d::Identity());
    Field<Mat3d> b(makeStandardBox(4), Mat3d::Identity());
    auto la = expr::field(a);
    auto lb = expr::field(b);
    try {
        (void)expr::doubleContract(la, lb);
        FAIL() << "Expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        const std::string what = e.what();
        EXPECT_NE(what.find("tensor::expr::DoubleContract"), std::string::npos);
    }
}

TEST(TensorGridMismatchTest, ErrorMessageNamesEntryPoint) {
    Field<Mat3d> a(makeStandardBox(8), Mat3d::Identity());
    Field<Mat3d> b(makeStandardBox(4), Mat3d::Identity());
    try {
        (void)tensor::doubleContract(a, b);
        FAIL() << "Expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        const std::string what = e.what();
        EXPECT_NE(what.find("tensor::doubleContract"), std::string::npos);
        EXPECT_NE(what.find("Nx"), std::string::npos);
    }
}

}  // namespace
