#include <gtest/gtest.h>

#include <string_view>

#include <gridcalc/version.hpp>

TEST(VersionTest, ReturnsZeroOneZero) {
    EXPECT_EQ(std::string_view{gridcalc::version()}, std::string_view{"0.1.0"});
}
