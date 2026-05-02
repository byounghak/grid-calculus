# Dependency pins for gridcalc.
#
# Per specs/tech-stack.md ("Dependency policy"):
#   - Header-only or vendor-able only.
#   - Pulled via CMake FetchContent with pinned tags / commit hashes.
#   - Updates are explicit changes reviewed in a PR.
#
# Resolved commit hashes are recorded next to each tag so a tag retag (rare but
# possible) cannot silently change the source tree under us.

include(FetchContent)

# -- Eigen ---------------------------------------------------------------------
# Tag: 3.4.0  (resolved commit: 3147391d946bb4b6c68edd901f2add6ac1f31f8c)
FetchContent_Declare(
    Eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3147391d946bb4b6c68edd901f2add6ac1f31f8c
    GIT_SHALLOW    TRUE
)
set(EIGEN_BUILD_DOC      OFF CACHE INTERNAL "")
set(EIGEN_BUILD_TESTING  OFF CACHE INTERNAL "")
set(BUILD_TESTING        OFF CACHE INTERNAL "")
set(EIGEN_BUILD_PKGCONFIG OFF CACHE INTERNAL "")

# -- GoogleTest ----------------------------------------------------------------
# Tag: v1.14.0  (resolved commit: f8d7d77c06936315286eb55f8de22cd23c188571)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        f8d7d77c06936315286eb55f8de22cd23c188571
    GIT_SHALLOW    TRUE
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(Eigen googletest)
