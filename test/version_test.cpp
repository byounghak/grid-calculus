#include <gtest/gtest.h>

#include <string_view>

#include <gridcalc/version.hpp>

TEST(VersionTest, ReturnsZeroFourteenFive) {
    EXPECT_EQ(std::string_view{gridcalc::version()}, std::string_view{"0.14.5"});
}
