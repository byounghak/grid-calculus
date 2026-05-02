include(GoogleTest)

add_executable(gridcalc_tests
    version_test.cpp
    index_policy_test.cpp
    grid_test.cpp
    field_test.cpp
    laplacian_test.cpp
)

target_link_libraries(gridcalc_tests PRIVATE
    gridcalc::gridcalc
    gridcalc_warnings
    GTest::gtest_main
)

gtest_discover_tests(gridcalc_tests)
