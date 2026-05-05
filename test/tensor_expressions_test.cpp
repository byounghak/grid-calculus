// Phase 15 — contraction expression-template (ET) layer.
//
// Each test in this file pairs a `tensor::expr::*` lazy node with the
// corresponding eager Phase 13 primitive (or a hand-rolled reference for
// the new `Sym` / `IdentityField` cases) and asserts byte-for-byte
// agreement after `materialize`. The fused `func::integrate(expr, tag)`
// overload is verified against the materialize-then-integrate path.

#include <cmath>

#include <gtest/gtest.h>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/func/Integrate.hpp>
#include <gridcalc/tensor/Contraction.hpp>
#include <gridcalc/tensor/Expressions.hpp>

namespace {

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Mat3d;
using gridcalc::core::Vec3d;
namespace expr = gridcalc::tensor::expr;
namespace tensor = gridcalc::tensor;

constexpr double kPi = 3.14159265358979323846;

Grid makeBox(int N) {
    const double h = 2.0 * kPi / N;
    return Grid(N, N, N, Vec3d(h, h, h));
}

// Sample a Mat3d-valued callable into a Field<Mat3d>.
template <class F>
Field<Mat3d> sampleMat(const Grid& g, F&& f) {
    Field<Mat3d> M(g, Mat3d::Zero());
    const auto& h = g.getCellSize();
    for (int k = 0; k < g.getNz(); ++k) {
        for (int j = 0; j < g.getNy(); ++j) {
            for (int i = 0; i < g.getNx(); ++i) {
                M(i, j, k) = f(static_cast<double>(i) * h(0),
                               static_cast<double>(j) * h(1),
                               static_cast<double>(k) * h(2));
            }
        }
    }
    return M;
}

// A non-trivial periodic Mat3d field with all 9 components active.
auto kManufacturedMatA = [](double x, double y, double z) {
    Mat3d m;
    m << std::sin(x), std::cos(y), std::sin(z),
         std::cos(x) * std::sin(y), std::sin(2.0 * x), std::cos(z),
         std::sin(x) * std::cos(z), std::cos(x) * std::sin(z), std::sin(x + y + z);
    return m;
};

auto kManufacturedMatB = [](double x, double y, double z) {
    Mat3d m;
    m << std::cos(x), std::sin(2.0 * y), std::cos(z + x),
         std::sin(x) * std::cos(y), std::cos(2.0 * z), std::sin(y),
         std::cos(y) * std::sin(z), std::sin(x) * std::cos(y), std::cos(x + y);
    return m;
};

// -------------------------------------------------------------------------
//  Leaf and IdentityField
// -------------------------------------------------------------------------

TEST(TensorExpressionsTest, LeafEvaluatesToFieldValue) {
    const Grid g = makeBox(8);
    const Field<Mat3d> M = sampleMat(g, kManufacturedMatA);
    const auto leaf = expr::field(M);
    for (int k = 0; k < g.getNz(); ++k) {
        for (int j = 0; j < g.getNy(); ++j) {
            for (int i = 0; i < g.getNx(); ++i) {
                EXPECT_EQ(leaf.evalAt(i, j, k), M(i, j, k));
            }
        }
    }
}

TEST(TensorExpressionsTest, IdentityFieldEvaluatesToMatIdentity) {
    const Grid g = makeBox(4);
    const auto id = expr::identityField(g);
    const Mat3d expected = Mat3d::Identity();
    for (int k = 0; k < g.getNz(); ++k) {
        for (int j = 0; j < g.getNy(); ++j) {
            for (int i = 0; i < g.getNx(); ++i) {
                EXPECT_EQ(id.evalAt(i, j, k), expected);
            }
        }
    }
}

// -------------------------------------------------------------------------
//  Combinator nodes — Scaled, Plus, Minus, Negate
// -------------------------------------------------------------------------

TEST(TensorExpressionsTest, ScaledMatchesEagerScalarMultiply) {
    const Grid g = makeBox(8);
    const Field<Mat3d> M = sampleMat(g, kManufacturedMatA);
    const double c = 2.5;
    const Field<Mat3d> result = expr::materialize(c * expr::field(M));
    for (int k = 0; k < g.getNz(); ++k) {
        for (int j = 0; j < g.getNy(); ++j) {
            for (int i = 0; i < g.getNx(); ++i) {
                const Mat3d expected = c * M(i, j, k);
                EXPECT_TRUE(result(i, j, k).isApprox(expected, 1e-15));
            }
        }
    }
}

TEST(TensorExpressionsTest, ScalarRightMultiplyAndDivisionMatchEager) {
    const Grid g = makeBox(8);
    const Field<Mat3d> M = sampleMat(g, kManufacturedMatA);
    const double c = 4.0;
    const Field<Mat3d> right = expr::materialize(expr::field(M) * c);
    const Field<Mat3d> div = expr::materialize(expr::field(M) / c);
    for (int k = 0; k < g.getNz(); ++k) {
        for (int j = 0; j < g.getNy(); ++j) {
            for (int i = 0; i < g.getNx(); ++i) {
                EXPECT_TRUE(right(i, j, k).isApprox(c * M(i, j, k), 1e-15));
                EXPECT_TRUE(div(i, j, k).isApprox(M(i, j, k) / c, 1e-15));
            }
        }
    }
}

TEST(TensorExpressionsTest, PlusMinusMatchesEagerSumDiff) {
    const Grid g = makeBox(6);
    const Field<Mat3d> A = sampleMat(g, kManufacturedMatA);
    const Field<Mat3d> B = sampleMat(g, kManufacturedMatB);
    const Field<Mat3d> sum = expr::materialize(expr::field(A) + expr::field(B));
    const Field<Mat3d> diff = expr::materialize(expr::field(A) - expr::field(B));
    for (int k = 0; k < g.getNz(); ++k) {
        for (int j = 0; j < g.getNy(); ++j) {
            for (int i = 0; i < g.getNx(); ++i) {
                EXPECT_TRUE(sum(i, j, k).isApprox(A(i, j, k) + B(i, j, k), 1e-15));
                EXPECT_TRUE(diff(i, j, k).isApprox(A(i, j, k) - B(i, j, k), 1e-15));
            }
        }
    }
}

TEST(TensorExpressionsTest, NegateMatchesEager) {
    const Grid g = makeBox(6);
    const Field<Mat3d> A = sampleMat(g, kManufacturedMatA);
    const Field<Mat3d> result = expr::materialize(-expr::field(A));
    for (int k = 0; k < g.getNz(); ++k) {
        for (int j = 0; j < g.getNy(); ++j) {
            for (int i = 0; i < g.getNx(); ++i) {
                EXPECT_TRUE(result(i, j, k).isApprox(-A(i, j, k), 1e-15));
            }
        }
    }
}

// -------------------------------------------------------------------------
//  Tensor primitives — Trace, Sym, SingleContract, DoubleContract
// -------------------------------------------------------------------------

TEST(TensorExpressionsTest, TraceMatchesEagerPhase13) {
    const Grid g = makeBox(8);
    const Field<Mat3d> A = sampleMat(g, kManufacturedMatA);
    const Field<double> lazy = expr::materialize(expr::trace(expr::field(A)));
    const Field<double> eager = tensor::trace(A);
    for (int k = 0; k < g.getNz(); ++k) {
        for (int j = 0; j < g.getNy(); ++j) {
            for (int i = 0; i < g.getNx(); ++i) {
                EXPECT_DOUBLE_EQ(lazy(i, j, k), eager(i, j, k));
            }
        }
    }
}

TEST(TensorExpressionsTest, SymMatchesHandComputed) {
    const Grid g = makeBox(6);
    const Field<Mat3d> A = sampleMat(g, kManufacturedMatA);
    const Field<Mat3d> result = expr::materialize(expr::sym(expr::field(A)));
    for (int k = 0; k < g.getNz(); ++k) {
        for (int j = 0; j < g.getNy(); ++j) {
            for (int i = 0; i < g.getNx(); ++i) {
                const Mat3d m = A(i, j, k);
                const Mat3d expected = 0.5 * (m + m.transpose());
                EXPECT_TRUE(result(i, j, k).isApprox(expected, 1e-15));
            }
        }
    }
}

TEST(TensorExpressionsTest, SingleContractMatchesEagerPhase13) {
    const Grid g = makeBox(6);
    const Field<Mat3d> A = sampleMat(g, kManufacturedMatA);
    const Field<Mat3d> B = sampleMat(g, kManufacturedMatB);
    const Field<Mat3d> lazy = expr::materialize(expr::singleContract(expr::field(A), expr::field(B)));
    const Field<Mat3d> eager = tensor::singleContract(A, B);
    for (int k = 0; k < g.getNz(); ++k) {
        for (int j = 0; j < g.getNy(); ++j) {
            for (int i = 0; i < g.getNx(); ++i) {
                EXPECT_TRUE(lazy(i, j, k).isApprox(eager(i, j, k), 1e-15));
            }
        }
    }
}

TEST(TensorExpressionsTest, DoubleContractMatchesEagerPhase13) {
    const Grid g = makeBox(6);
    const Field<Mat3d> A = sampleMat(g, kManufacturedMatA);
    const Field<Mat3d> B = sampleMat(g, kManufacturedMatB);
    const Field<double> lazy = expr::materialize(expr::doubleContract(expr::field(A), expr::field(B)));
    const Field<double> eager = tensor::doubleContract(A, B);
    for (int k = 0; k < g.getNz(); ++k) {
        for (int j = 0; j < g.getNy(); ++j) {
            for (int i = 0; i < g.getNx(); ++i) {
                EXPECT_DOUBLE_EQ(lazy(i, j, k), eager(i, j, k));
            }
        }
    }
}

// -------------------------------------------------------------------------
//  Composition tests
// -------------------------------------------------------------------------

TEST(TensorExpressionsTest, NestedExpression_TraceOfSym) {
    const Grid g = makeBox(8);
    const Field<Mat3d> A = sampleMat(g, kManufacturedMatA);
    const Field<double> result = expr::materialize(expr::trace(expr::sym(expr::field(A))));
    for (int k = 0; k < g.getNz(); ++k) {
        for (int j = 0; j < g.getNy(); ++j) {
            for (int i = 0; i < g.getNx(); ++i) {
                const Mat3d m = A(i, j, k);
                const Mat3d s = 0.5 * (m + m.transpose());
                EXPECT_DOUBLE_EQ(result(i, j, k), s.trace());
            }
        }
    }
}

TEST(TensorExpressionsTest, LinearCombinationOfDoubleContracts) {
    const Grid g = makeBox(6);
    const Field<Mat3d> A = sampleMat(g, kManufacturedMatA);
    const Field<Mat3d> B = sampleMat(g, kManufacturedMatB);
    const double c1 = 0.5;
    const double c2 = 1.25;
    auto leafA = expr::field(A);
    auto leafB = expr::field(B);
    const Field<double> lazy = expr::materialize(
        c1 * expr::doubleContract(leafA, leafB) + c2 * expr::doubleContract(leafB, leafA));
    for (int k = 0; k < g.getNz(); ++k) {
        for (int j = 0; j < g.getNy(); ++j) {
            for (int i = 0; i < g.getNx(); ++i) {
                const Mat3d a = A(i, j, k);
                const Mat3d b = B(i, j, k);
                const double ab = (a.array() * b.array()).sum();
                const double ba = (b.array() * a.array()).sum();
                const double expected = c1 * ab + c2 * ba;
                EXPECT_DOUBLE_EQ(lazy(i, j, k), expected);
            }
        }
    }
}

TEST(TensorExpressionsTest, IdentityFieldComposesWithDoubleContract) {
    // doubleContract(I, M) = trace(M) at every cell, since I_ij M_ij = M_ii.
    // The two paths sum the same three diagonal entries but in different
    // orders (the doubleContract path adds six explicit zeros), so the
    // identity holds only up to round-off, not bit-exactly.
    const Grid g = makeBox(6);
    const Field<Mat3d> A = sampleMat(g, kManufacturedMatA);
    const Field<double> via_id = expr::materialize(
        expr::doubleContract(expr::identityField(g), expr::field(A)));
    const Field<double> via_trace = tensor::trace(A);
    for (int k = 0; k < g.getNz(); ++k) {
        for (int j = 0; j < g.getNy(); ++j) {
            for (int i = 0; i < g.getNx(); ++i) {
                EXPECT_NEAR(via_id(i, j, k), via_trace(i, j, k), 1e-14);
            }
        }
    }
}

TEST(TensorExpressionsTest, IdentityFieldComposesWithSingleContract) {
    // singleContract(I, M) = M at every cell.
    const Grid g = makeBox(6);
    const Field<Mat3d> A = sampleMat(g, kManufacturedMatA);
    const Field<Mat3d> result = expr::materialize(
        expr::singleContract(expr::identityField(g), expr::field(A)));
    for (int k = 0; k < g.getNz(); ++k) {
        for (int j = 0; j < g.getNy(); ++j) {
            for (int i = 0; i < g.getNx(); ++i) {
                EXPECT_TRUE(result(i, j, k).isApprox(A(i, j, k), 1e-15));
            }
        }
    }
}

// -------------------------------------------------------------------------
//  Fused integrate
// -------------------------------------------------------------------------

TEST(TensorExpressionsTest, IntegrateExpr_FusedReduction) {
    const Grid g = makeBox(8);
    const Field<Mat3d> A = sampleMat(g, kManufacturedMatA);
    const Field<Mat3d> B = sampleMat(g, kManufacturedMatB);
    auto leafA = expr::field(A);
    auto leafB = expr::field(B);
    const double via_expr =
        gridcalc::func::integrate(0.5 * expr::doubleContract(leafA, leafB));
    const Field<double> integrand =
        expr::materialize(0.5 * expr::doubleContract(leafA, leafB));
    const double via_field = gridcalc::func::integrate(integrand);
    EXPECT_DOUBLE_EQ(via_expr, via_field);
}

TEST(TensorExpressionsTest, IntegrateExpr_KahanReduction) {
    const Grid g = makeBox(8);
    const Field<Mat3d> A = sampleMat(g, kManufacturedMatA);
    const Field<Mat3d> B = sampleMat(g, kManufacturedMatB);
    auto leafA = expr::field(A);
    auto leafB = expr::field(B);
    const double via_expr = gridcalc::func::integrate(
        expr::doubleContract(leafA, leafB), gridcalc::func::Kahan{});
    const Field<double> integrand =
        expr::materialize(expr::doubleContract(leafA, leafB));
    const double via_field =
        gridcalc::func::integrate(integrand, gridcalc::func::Kahan{});
    EXPECT_DOUBLE_EQ(via_expr, via_field);
}

}  // namespace
