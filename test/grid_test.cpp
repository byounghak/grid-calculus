#include <gtest/gtest.h>

#include <stdexcept>

#include <gridcalc/core/Grid.hpp>

using gridcalc::core::Grid;
using gridcalc::core::Vec3d;

namespace {

constexpr double kPi = 3.14159265358979323846;

}  // namespace

TEST(GridTest, ConstructsAndExposesShape) {
  Grid g(8, 16, 32, Vec3d(0.5, 0.25, 0.125));
  EXPECT_EQ(g.getNx(), 8);
  EXPECT_EQ(g.getNy(), 16);
  EXPECT_EQ(g.getNz(), 32);
  EXPECT_DOUBLE_EQ(g.getCellSize()(0), 0.5);
  EXPECT_DOUBLE_EQ(g.getCellSize()(1), 0.25);
  EXPECT_DOUBLE_EQ(g.getCellSize()(2), 0.125);
  EXPECT_EQ(g.getSize(), 8u * 16u * 32u);
  EXPECT_EQ(g.size(), g.getSize());
}

TEST(GridTest, ComputesCellVolume) {
  Grid g(4, 4, 4, Vec3d(2.0, 3.0, 5.0));
  EXPECT_DOUBLE_EQ(g.getCellVolume(), 30.0);
}

TEST(GridTest, LinearIndexFormulaIsIFastest) {
  const int Nx = 4;
  const int Ny = 5;
  const int Nz = 6;
  Grid g(Nx, Ny, Nz, Vec3d(1.0, 1.0, 1.0));

  EXPECT_EQ(g.getLinearIndex(0, 0, 0), 0u);
  EXPECT_EQ(g.getLinearIndex(1, 0, 0), 1u);
  EXPECT_EQ(g.getLinearIndex(Nx - 1, 0, 0), static_cast<std::size_t>(Nx - 1));
  EXPECT_EQ(g.getLinearIndex(0, 1, 0), static_cast<std::size_t>(Nx));
  EXPECT_EQ(g.getLinearIndex(0, 0, 1), static_cast<std::size_t>(Nx) * Ny);
  EXPECT_EQ(g.getLinearIndex(2, 3, 4),
            static_cast<std::size_t>(2 + Nx * (3 + Ny * 4)));
}

TEST(GridTest, LinearIndexCoversWholeStorage) {
  const int Nx = 3;
  const int Ny = 4;
  const int Nz = 5;
  Grid g(Nx, Ny, Nz, Vec3d(1.0, 1.0, 1.0));

  std::size_t expected = 0;
  for (int k = 0; k < Nz; ++k) {
    for (int j = 0; j < Ny; ++j) {
      for (int i = 0; i < Nx; ++i) {
        EXPECT_EQ(g.getLinearIndex(i, j, k), expected);
        ++expected;
      }
    }
  }
  EXPECT_EQ(expected, g.getSize());
}

TEST(GridTest, RejectsZeroOrNegativeShape) {
  EXPECT_THROW(Grid(0, 4, 4, Vec3d(1, 1, 1)), std::invalid_argument);
  EXPECT_THROW(Grid(4, 0, 4, Vec3d(1, 1, 1)), std::invalid_argument);
  EXPECT_THROW(Grid(4, 4, 0, Vec3d(1, 1, 1)), std::invalid_argument);
  EXPECT_THROW(Grid(-1, 4, 4, Vec3d(1, 1, 1)), std::invalid_argument);
}

TEST(GridTest, RejectsNonPositiveSpacing) {
  EXPECT_THROW(Grid(4, 4, 4, Vec3d(0.0, 1.0, 1.0)), std::invalid_argument);
  EXPECT_THROW(Grid(4, 4, 4, Vec3d(1.0, -1.0, 1.0)), std::invalid_argument);
  EXPECT_THROW(Grid(4, 4, 4, Vec3d(1.0, 1.0, 0.0)), std::invalid_argument);
}

TEST(GridTest, AcceptsNonCubicAndAnisotropicSpacing) {
  EXPECT_NO_THROW(Grid(4, 8, 16, Vec3d(2.0 * kPi / 4, 2.0 * kPi / 8, 2.0 * kPi / 16)));
}
