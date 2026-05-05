// Tests for `gridcalc::core::Grid::operator==` (added at 0.14.2).
//
// Bit-exact componentwise comparison of (Nx, Ny, Nz) and the three
// components of cell_size. Used by the time-integration machinery
// (solve::integrate, detail::axpyFresh) to validate that an RHS
// callable returns a Field allocated against the same Grid as the
// integrator's state.

#include <gtest/gtest.h>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Grid.hpp>

using gridcalc::core::Grid;
using gridcalc::core::Vec3d;

TEST(GridEqualityTest, IdenticalGridsCompareEqual) {
  Grid a(8, 8, 8, Vec3d(0.1, 0.2, 0.3));
  Grid b(8, 8, 8, Vec3d(0.1, 0.2, 0.3));
  EXPECT_TRUE(a == b);
  EXPECT_FALSE(a != b);
}

TEST(GridEqualityTest, ReflexiveAndCopy) {
  Grid a(8, 8, 8, Vec3d(0.1, 0.2, 0.3));
  EXPECT_TRUE(a == a);
  Grid b = a;
  EXPECT_TRUE(a == b);
  EXPECT_TRUE(b == a);
}

TEST(GridEqualityTest, DiffersOnNx) {
  Grid a(8, 8, 8, Vec3d(0.1, 0.1, 0.1));
  Grid b(7, 8, 8, Vec3d(0.1, 0.1, 0.1));
  EXPECT_FALSE(a == b);
  EXPECT_TRUE(a != b);
}

TEST(GridEqualityTest, DiffersOnNy) {
  Grid a(8, 8, 8, Vec3d(0.1, 0.1, 0.1));
  Grid b(8, 9, 8, Vec3d(0.1, 0.1, 0.1));
  EXPECT_FALSE(a == b);
}

TEST(GridEqualityTest, DiffersOnNz) {
  Grid a(8, 8, 8, Vec3d(0.1, 0.1, 0.1));
  Grid b(8, 8, 16, Vec3d(0.1, 0.1, 0.1));
  EXPECT_FALSE(a == b);
}

TEST(GridEqualityTest, DiffersOnHx) {
  Grid a(8, 8, 8, Vec3d(0.1, 0.1, 0.1));
  Grid b(8, 8, 8, Vec3d(0.2, 0.1, 0.1));
  EXPECT_FALSE(a == b);
}

TEST(GridEqualityTest, DiffersOnHy) {
  Grid a(8, 8, 8, Vec3d(0.1, 0.1, 0.1));
  Grid b(8, 8, 8, Vec3d(0.1, 0.2, 0.1));
  EXPECT_FALSE(a == b);
}

TEST(GridEqualityTest, DiffersOnHz) {
  Grid a(8, 8, 8, Vec3d(0.1, 0.1, 0.1));
  Grid b(8, 8, 8, Vec3d(0.1, 0.1, 0.2));
  EXPECT_FALSE(a == b);
}

TEST(GridEqualityTest, BitExactNotTolerance) {
  // Document the bit-exact contract: 0.1 + 0.2 != 0.3 at the bit level
  // (IEEE 754 rounding). Two grids constructed with arithmetically-
  // equal but bit-distinct cell sizes do NOT compare equal. This is by
  // design -- there is no source of numerical drift on Grid values in
  // this codebase, so any non-bit-exact mismatch is by definition a
  // different Grid object.
  const double drifted = 0.1 + 0.2;  // 0.30000000000000004 in IEEE 754
  ASSERT_NE(drifted, 0.3);  // sanity-check the platform double rounding
  Grid a(8, 8, 8, Vec3d(0.3, 0.1, 0.1));
  Grid b(8, 8, 8, Vec3d(drifted, 0.1, 0.1));
  EXPECT_FALSE(a == b)
      << "bit-exact contract: arithmetically-equal but bit-distinct "
         "cell sizes must compare unequal";
}
