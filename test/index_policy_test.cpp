#include <gtest/gtest.h>

#include <gridcalc/core/IndexPolicy.hpp>

using gridcalc::core::IndexPolicy::Periodic;

TEST(PeriodicTest, WrapsInRange) {
  for (int i = 0; i < 5; ++i) {
    EXPECT_EQ(Periodic::wrap(i, 5), i);
  }
}

TEST(PeriodicTest, WrapsNegative) {
  EXPECT_EQ(Periodic::wrap(-1, 5), 4);
  EXPECT_EQ(Periodic::wrap(-5, 5), 0);
  EXPECT_EQ(Periodic::wrap(-7, 5), 3);
  EXPECT_EQ(Periodic::wrap(-13, 5), 2);
}

TEST(PeriodicTest, WrapsLarge) {
  EXPECT_EQ(Periodic::wrap(5, 5), 0);
  EXPECT_EQ(Periodic::wrap(6, 5), 1);
  EXPECT_EQ(Periodic::wrap(13, 5), 3);
  EXPECT_EQ(Periodic::wrap(100, 5), 0);
}

TEST(PeriodicTest, WrapsSingleCellExtent) {
  EXPECT_EQ(Periodic::wrap(0, 1), 0);
  EXPECT_EQ(Periodic::wrap(7, 1), 0);
  EXPECT_EQ(Periodic::wrap(-3, 1), 0);
}

TEST(PeriodicTest, IsConstexpr) {
  constexpr int wrapped = Periodic::wrap(-1, 4);
  static_assert(wrapped == 3, "Periodic::wrap must be usable in constexpr context");
  EXPECT_EQ(wrapped, 3);
}
