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
    spectral_basic_test.cpp
    fd_fft_crosscheck_test.cpp
    pfc_functional_test.cpp
    vtk_writer_test.cpp
    cahn_hilliard_test.cpp
    vector_tensor_test.cpp
)

# Tests for the example helpers (`examples/common/...`) include them
# directly as headers; they are not gated on GRIDCALC_BUILD_EXAMPLES so
# that the round-trip and acceptance gates run on every PR.
target_include_directories(gridcalc_tests PRIVATE
    "${CMAKE_SOURCE_DIR}/examples"
)

target_link_libraries(gridcalc_tests PRIVATE
    gridcalc::gridcalc
    gridcalc_warnings
    GTest::gtest_main
)

gtest_discover_tests(gridcalc_tests)
