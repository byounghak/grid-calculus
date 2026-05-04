#include <gtest/gtest.h>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Mat3d;
using gridcalc::core::Vec3d;

TEST(FieldTest, ZeroInitFillsZero) {
  Grid g(4, 5, 6, Vec3d(1.0, 1.0, 1.0));
  Field<double> f(g);
  EXPECT_EQ(f.getSize(), g.getSize());
  for (double v : f) {
    EXPECT_DOUBLE_EQ(v, 0.0);
  }
}

TEST(FieldTest, ValueInitBroadcasts) {
  Grid g(3, 3, 3, Vec3d(1.0, 1.0, 1.0));
  Field<double> f(g, 7.5);
  for (double v : f) {
    EXPECT_DOUBLE_EQ(v, 7.5);
  }
}

TEST(FieldTest, CallableInitSamplesCartesian) {
  Grid g(4, 4, 4, Vec3d(0.5, 1.0, 2.0));
  Field<double> f(g, [](double x, double y, double z) { return x + 2.0 * y + 3.0 * z; });

  for (int k = 0; k < g.getNz(); ++k) {
    for (int j = 0; j < g.getNy(); ++j) {
      for (int i = 0; i < g.getNx(); ++i) {
        const double expected =
            i * 0.5 + 2.0 * (j * 1.0) + 3.0 * (k * 2.0);
        EXPECT_DOUBLE_EQ(f(i, j, k), expected);
      }
    }
  }
}

TEST(FieldTest, PeriodicWrapAtNegativeIndex) {
  Grid g(4, 4, 4, Vec3d(1.0, 1.0, 1.0));
  Field<double> f(g, [](double x, double y, double z) { return 100.0 * x + 10.0 * y + z; });
  EXPECT_DOUBLE_EQ(f(-1, 0, 0), f(g.getNx() - 1, 0, 0));
  EXPECT_DOUBLE_EQ(f(0, -1, 0), f(0, g.getNy() - 1, 0));
  EXPECT_DOUBLE_EQ(f(0, 0, -1), f(0, 0, g.getNz() - 1));
  EXPECT_DOUBLE_EQ(f(-7, -3, -2), f(g.getNx() - 3, g.getNy() - 3, g.getNz() - 2));
}

TEST(FieldTest, PeriodicWrapAtOverflowIndex) {
  Grid g(4, 4, 4, Vec3d(1.0, 1.0, 1.0));
  Field<double> f(g, [](double x, double y, double z) { return 100.0 * x + 10.0 * y + z; });
  EXPECT_DOUBLE_EQ(f(g.getNx(), 0, 0), f(0, 0, 0));
  EXPECT_DOUBLE_EQ(f(0, g.getNy() + 1, 0), f(0, 1, 0));
  EXPECT_DOUBLE_EQ(f(0, 0, g.getNz() + 2), f(0, 0, 2));
}

TEST(FieldTest, WriteThenReadRoundTrip) {
  Grid g(2, 3, 4, Vec3d(1.0, 1.0, 1.0));
  Field<double> f(g);
  for (int k = 0; k < g.getNz(); ++k) {
    for (int j = 0; j < g.getNy(); ++j) {
      for (int i = 0; i < g.getNx(); ++i) {
        f(i, j, k) = static_cast<double>(100 * i + 10 * j + k);
      }
    }
  }
  for (int k = 0; k < g.getNz(); ++k) {
    for (int j = 0; j < g.getNy(); ++j) {
      for (int i = 0; i < g.getNx(); ++i) {
        EXPECT_DOUBLE_EQ(f(i, j, k), static_cast<double>(100 * i + 10 * j + k));
      }
    }
  }
}

TEST(FieldTest, GridCopiedByValue) {
  Grid g(2, 2, 2, Vec3d(1.0, 1.0, 1.0));
  Field<double> f(g);
  EXPECT_EQ(f.getGrid().getNx(), g.getNx());
  EXPECT_EQ(f.getGrid().getNy(), g.getNy());
  EXPECT_EQ(f.getGrid().getNz(), g.getNz());
  EXPECT_DOUBLE_EQ(f.getGrid().getCellVolume(), g.getCellVolume());
}

TEST(FieldMat3dTest, BroadcastConstructorFillsEveryCell) {
  Grid g(3, 4, 5, Vec3d(1.0, 1.0, 1.0));
  Mat3d M;
  M << 1.0, 2.0, 3.0,
       4.0, 5.0, 6.0,
       7.0, 8.0, 9.0;
  Field<Mat3d> f(g, M);
  EXPECT_EQ(f.getSize(), g.getSize());
  for (int k = 0; k < g.getNz(); ++k) {
    for (int j = 0; j < g.getNy(); ++j) {
      for (int i = 0; i < g.getNx(); ++i) {
        EXPECT_EQ(f(i, j, k), M);
      }
    }
  }
}

TEST(FieldMat3dTest, ZeroBroadcastIsZeroEverywhere) {
  Grid g(2, 3, 2, Vec3d(1.0, 1.0, 1.0));
  Field<Mat3d> f(g, Mat3d::Zero());
  for (auto& M : f) {
    EXPECT_TRUE(M.isZero());
  }
}

TEST(FieldMat3dTest, WriteThenReadRoundTrip) {
  Grid g(2, 3, 4, Vec3d(1.0, 1.0, 1.0));
  Field<Mat3d> f(g, Mat3d::Zero());
  for (int k = 0; k < g.getNz(); ++k) {
    for (int j = 0; j < g.getNy(); ++j) {
      for (int i = 0; i < g.getNx(); ++i) {
        Mat3d M = Mat3d::Zero();
        M(0, 0) = static_cast<double>(100 * i + 10 * j + k);
        M(1, 2) = static_cast<double>(i + j + k);
        f(i, j, k) = M;
      }
    }
  }
  for (int k = 0; k < g.getNz(); ++k) {
    for (int j = 0; j < g.getNy(); ++j) {
      for (int i = 0; i < g.getNx(); ++i) {
        EXPECT_DOUBLE_EQ(f(i, j, k)(0, 0), static_cast<double>(100 * i + 10 * j + k));
        EXPECT_DOUBLE_EQ(f(i, j, k)(1, 2), static_cast<double>(i + j + k));
      }
    }
  }
}

TEST(FieldMat3dTest, IFastestLayoutMatchesScalarField) {
  // Storage layout for Field<Mat3d> uses the same Grid::getLinearIndex
  // mapping as Field<double>, so adjacent grid points along x are
  // adjacent in `data()` (one Mat3d apart, not one double apart).
  Grid g(4, 2, 2, Vec3d(1.0, 1.0, 1.0));
  Field<Mat3d> f(g, Mat3d::Zero());
  for (int k = 0; k < g.getNz(); ++k) {
    for (int j = 0; j < g.getNy(); ++j) {
      for (int i = 0; i < g.getNx(); ++i) {
        f(i, j, k) = static_cast<double>(i) * Mat3d::Identity();
      }
    }
  }
  // i-fastest: (1, 0, 0) - (0, 0, 0) is one Mat3d step in storage.
  const Mat3d* base = f.data();
  EXPECT_EQ(base[0], 0.0 * Mat3d::Identity());
  EXPECT_EQ(base[1], 1.0 * Mat3d::Identity());
  EXPECT_EQ(base[2], 2.0 * Mat3d::Identity());
  EXPECT_EQ(base[3], 3.0 * Mat3d::Identity());
}
