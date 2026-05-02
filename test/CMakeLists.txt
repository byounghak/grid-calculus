include(GoogleTest)

add_executable(gridcalc_tests
    version_test.cpp
)

target_link_libraries(gridcalc_tests PRIVATE
    gridcalc::gridcalc
    gridcalc_warnings
    GTest::gtest_main
)

gtest_discover_tests(gridcalc_tests)
