# Validation: Fix silent stencil aliasing on small periodic axes

## Build

- [ ] CI green on Ubuntu GCC.
- [ ] CI green on Ubuntu Clang.
- [ ] CI green on `clang-debug-nofft` configuration (the precondition is unrelated to the FFT cross-check; the build must stay green either way).
- [ ] Warnings clean under the project's strict warning set (`-Wall -Wextra -Wpedantic -Wconversion`).
- [ ] `clang-format` and `clang-tidy` pass.
- [ ] `docs-build` job green (Doxygen + LaTeX User Guide + Developer Note).
- [ ] `docs-lint` job green (Doxygen `WARN_IF_UNDOCUMENTED` + `\since`-tag regex; the new helper carries `\since 0.14.1` and the touched operators carry an updated `\throws` line on their existing Doxygen blocks).

## Tests

- [ ] `test/stencil_axis_extent_test.cpp` exists and asserts the throw at `N = 2 * radius` for each of `radius ∈ {1, 2, 3}`.
- [ ] Each axis (x, y, z) is independently exercised to confirm the helper reports the correct axis label in its error message.
- [ ] At least one test confirms that an axis with per-axis derivative count `0` does **not** trip the precondition even at small extent, through `diff::detail::mixedDerivative`.
- [ ] Degenerate-axis carve-out tests: `N = 1` is accepted across `laplacian<4>`, `gradient<4>`, `biharmonic<4>`, `fvm::cellLaplacian` (the four representative operators spanning radius 2, 2, 3, 1 respectively). Each asserts no-throw and finite output.
- [ ] Smallest-valid-N boundary tests at `N = 2 * radius + 1` for every Phase 7 / Phase 8 / Phase 13 / Phase 14 operator (Laplacian, scalar / vector gradient, scalar / tensor divergence, biharmonic, mixed-partial Order=2 and Order=4, FVM `cellLaplacian`), each asserting no-throw and finite output.
- [ ] `test/fd_fft_crosscheck_test.cpp` is unmodified — the existing `{16, 32, 64, 128}` sweep stays the slope-analysis gate; the new file owns the boundary cases.
- [ ] Full test suite passes on `clang-debug` (existing 258 tests + new tests, all green).
- [ ] Full test suite passes on `clang-debug-nofft`.

## Documentation

- [ ] `\throws std::invalid_argument` line added to the Doxygen blocks of every touched operator (Laplacian, Gradient + Phase 13 overload, Divergence + Phase 13 overload, Biharmonic, every Phase 8 d-prefix partial via the shared `mixedDerivative` Doxygen, `fvm::cellLaplacian`, `fvm::diffuse`'s heterogeneous-D overload).
- [ ] New helper `gridcalc::diff::detail::requireAxisExtent` carries a Doxygen block with `\brief`, `\param`, `\throws`, `\since 0.14.1`. Lives under `detail/` so it is excluded from the public-API `\since` lint regex.
- [ ] `CHANGELOG.md` carries a new `0.14.1` block above the existing `0.14.0` block with `### Fixed`, `### Tests`, `### Internal` sections.
- [ ] `specs/STATUS.md` "Last updated" line bumped; "Phase 14 follow-up: PR #19" entry added describing this fix.

## Versioning

- [ ] `CMakeLists.txt` project `VERSION` is `0.14.1`.

## Manual verification

- [ ] On a fresh build, calling `diff::laplacian<4>(field)` with a grid where any axis has `N ≤ 4` raises `std::invalid_argument` with a message naming the axis, its `N`, the radius, and the required minimum.
- [ ] The same call on the same field with `Nx = Ny = Nz = 5` (the smallest valid case) succeeds and produces a numerically correct result (covered by the widened FD–FFT cross-check sweep).
- [ ] On uniform Cartesian grids with `N ≥ 5`, `diff::laplacian<4>` and `fvm::cellLaplacian` produce byte-for-byte the same output they did at `0.14.0` (no behavior change inside the validity domain).
