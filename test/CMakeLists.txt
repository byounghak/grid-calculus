include(GoogleTest)

add_executable(gridcalc_tests
    version_test.cpp
    index_policy_test.cpp
    grid_test.cpp
    field_test.cpp
    laplacian_test.cpp
    gradient_test.cpp
    divergence_test.cpp
    integrate_test.cpp
    functional_test.cpp
    diffusion_test.cpp
    integrator_test.cpp
    higher_order_test.cpp
    phase_8_test.cpp
)

target_link_libraries(gridcalc_tests PRIVATE
    gridcalc::gridcalc
    gridcalc_warnings
    GTest::gtest_main
)

gtest_discover_tests(gridcalc_tests)
