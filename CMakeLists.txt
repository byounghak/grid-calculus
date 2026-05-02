cmake_minimum_required(VERSION 3.20)

project(gridcalc
    VERSION 0.5.0
    DESCRIPTION "Differential operators and functional evaluation on periodic crystalline grids"
    LANGUAGES CXX
)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR
        "In-source builds are not allowed. "
        "Create a build directory and run cmake from there, e.g.:\n"
        "  cmake -B build && cmake --build build\n"
        "or use one of the presets: cmake --preset gcc-debug")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(GRIDCALC_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(GRIDCALC_BUILD_TESTS "Build the gridcalc test suite" ON)

add_library(gridcalc_warnings INTERFACE)
if(MSVC)
    target_compile_options(gridcalc_warnings INTERFACE /W4 /permissive-)
    if(GRIDCALC_WARNINGS_AS_ERRORS)
        target_compile_options(gridcalc_warnings INTERFACE /WX)
    endif()
else()
    target_compile_options(gridcalc_warnings INTERFACE
        -Wall -Wextra -Wpedantic -Wconversion)
    if(GRIDCALC_WARNINGS_AS_ERRORS)
        target_compile_options(gridcalc_warnings INTERFACE -Werror)
    endif()
endif()

include(cmake/Dependencies.cmake)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/include/gridcalc/version.hpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/generated/include/gridcalc/version.hpp"
    @ONLY
)

add_library(gridcalc INTERFACE)
add_library(gridcalc::gridcalc ALIAS gridcalc)
target_include_directories(gridcalc INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated/include>"
)
target_compile_features(gridcalc INTERFACE cxx_std_17)

# Eigen's own headers do not compile clean under -Wconversion. Inject its
# include path with the SYSTEM keyword so warnings from <Eigen/...> are
# suppressed while our own code still gets the full set of warnings.
get_target_property(_eigen_includes Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(gridcalc SYSTEM INTERFACE ${_eigen_includes})

if(GRIDCALC_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()
