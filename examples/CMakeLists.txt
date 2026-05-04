# Example programs that exercise the public API end-to-end.
#
# Targets here are built only when GRIDCALC_BUILD_EXAMPLES=ON. They are not
# part of the library install, are not run by CI (CI only type-checks
# them), and live alongside the library so they can call the same headers
# as production users.

add_executable(cahn_hilliard cahn_hilliard.cpp)

target_include_directories(cahn_hilliard PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(cahn_hilliard PRIVATE
    gridcalc::gridcalc
    gridcalc_warnings
)

add_executable(heterogeneous_diffusion heterogeneous_diffusion.cpp)

target_include_directories(heterogeneous_diffusion PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(heterogeneous_diffusion PRIVATE
    gridcalc::gridcalc
    gridcalc_warnings
)
