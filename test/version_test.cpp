#include <gtest/gtest.h>

#include <string_view>

#include <gridcalc/version.hpp>

TEST(VersionTest, ReturnsZeroTwoZero) {
    EXPECT_EQ(std::string_view{gridcalc::version()}, std::string_view{"0.2.0"});
}
