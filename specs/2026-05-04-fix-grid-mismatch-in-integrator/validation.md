# Validation: Fix silent Grid-mismatch on RHS return in time integrators

## Build

- [ ] CI green on Ubuntu GCC.
- [ ] CI green on Ubuntu Clang.
- [ ] CI green on the `clang-debug-nofft` configuration.
- [ ] Warnings clean under the project's strict warning set (`-Wall -Wextra -Wpedantic -Wconversion`).
- [ ] `clang-format` and `clang-tidy` pass.
- [ ] `docs-build` job green (Doxygen + LaTeX User Guide + Developer Note).
- [ ] `docs-lint` job green. New public symbol `Grid::operator==` carries `\since 0.14.2`; `solve::detail::requireSameGrid` lives under `detail/` and is excluded from the public-API `\since` regex.

## Tests

- [ ] `test/grid_equality_test.cpp` covers identical-grid equality, per-axis extent mismatch, per-component cell-size mismatch, and copy-equality. Bit-exact comparison verified (no tolerance papering over double mismatches).
- [ ] `test/integrator_grid_mismatch_test.cpp` covers:
  - [ ] `integrate(... RK4)` throws on a mismatched-RHS callable with different total cell count.
  - [ ] `integrate(... RK4)` throws on a mismatched-RHS callable with same total cell count but different per-axis extents (silent-layout-drift case).
  - [ ] `integrate(... RK4)` throws on a mismatched-RHS callable with same shape but different cell size (silent-`1/h`-drift case).
  - [ ] `integrate(... ExplicitEuler)` throws on each of the three mismatch flavors above.
  - [ ] `solve::detail::axpyFresh` throws when its two field arguments carry different grids.
  - [ ] Per-stage label sweep on RK4: the error message names `k1` / `k2` / `k3` / `k4` correctly, depending on which stage's RHS returned the bad grid.
  - [ ] Smoke test: a correctly-grid-preserving RHS (e.g., `[](const Field<double>& f){ return diff::laplacian(f); }`) integrates without throwing through both `RK4` and `ExplicitEuler`.
- [ ] Full test suite passes on `clang-debug` (existing 292 tests + new tests, all green).
- [ ] Full test suite passes on `clang-debug-nofft`.

## Documentation

- [ ] `\throws std::invalid_argument` line added to the Doxygen blocks of `integrate(... RK4)`, `integrate(... ExplicitEuler)`, `explicitEuler` (the Phase 5 wrapper), and `detail::axpyFresh`. Each names the new "RHS-grid mismatch" failure mode.
- [ ] `Grid::operator==` and `Grid::operator!=` carry Doxygen with `\brief`, `\returns`, `\since 0.14.2`.
- [ ] `solve::detail::requireSameGrid` carries Doxygen with `\brief`, `\param`, `\throws`, `\since 0.14.2`.
- [ ] `CHANGELOG.md` carries a new `0.14.2` block above the existing `0.14.1` block with `### Fixed`, `### Tests`, `### Internal` sections.
- [ ] `specs/STATUS.md` "Last updated" line bumped; "Phase 14 follow-up: PR #20" entry added describing this fix; the prior PR #19 entry is retained.

## Versioning

- [ ] `CMakeLists.txt` project `VERSION` is `0.14.2`.
- [ ] `test/version_test.cpp` asserts `"0.14.2"`.

## Manual verification

- [ ] An RHS callable that captures a different `Grid` and returns a `Field` allocated on it raises `std::invalid_argument` from `solve::integrate(...)` with a message naming the offending stage and the expected vs. actual `(Nx, Ny, Nz, hx, hy, hz)`.
- [ ] The existing FVM diffusion / Cahn–Hilliard / heterogeneous-D acceptance tests continue to pass without modification (the precondition does not regress the happy path).
- [ ] A `Grid` constructed with `Vec3d{0.1, 0.1, 0.1}` does **not** compare equal to one constructed with `Vec3d{0.05 + 0.05, 0.1, 0.1}` if those values differ at the bit level (documents the bit-exact contract for future readers).
