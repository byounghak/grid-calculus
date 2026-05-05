# Plan: Fix silent Grid-mismatch on RHS return in time integrators

## Group 1 — `Grid::operator==` / `operator!=`

- Add `bool operator==(const Grid&) const noexcept` and `operator!=` to [include/gridcalc/core/Grid.hpp](../../include/gridcalc/core/Grid.hpp). Bit-exact componentwise comparison of `_Nx, _Ny, _Nz, _cell_size`. Doxygen with `\since 0.14.2`.
- New `test/grid_equality_test.cpp` covering: identical grids compare equal; differing extent (each axis) compares unequal; differing cell size (each component) compares unequal; copy of a grid compares equal; bit-exact behavior — `0.1 + 0.2 != 0.3` style mismatch is not papered over (no tolerance applied).

**Commit message:** `feat(core): add Grid::operator== for bit-exact grid identity checks`

## Group 2 — precondition helper + apply to `axpyFresh` + integrators

- New header `include/gridcalc/solve/detail/RhsGridCheck.hpp` with the free function `requireSameGrid(const Field<double>& a, const Field<double>& b, const char* context)` that throws `std::invalid_argument` with a precise message (context + expected `(Nx, Ny, Nz, hx, hy, hz)` + actual) when `a.getGrid() != b.getGrid()`. Doxygen with `\since 0.14.2`.
- Apply at:
  - `detail::axpyFresh` in [include/gridcalc/solve/RK4.hpp](../../include/gridcalc/solve/RK4.hpp) — check `a.getGrid() == b.getGrid()` at entry.
  - `integrate(... RK4)` per stage — check that `k1`, `k2`, `k3`, `k4` each carry `psi.getGrid()` immediately after the corresponding `rhs(...)` call. Context strings name the stage.
  - `integrate(... ExplicitEuler)` in [include/gridcalc/solve/ExplicitEuler.hpp](../../include/gridcalc/solve/ExplicitEuler.hpp) — check `tmp.getGrid() == psi.getGrid()` after each `rhs(psi)` call.
- Update each touched function's Doxygen `\throws` line to reference the new precondition.

**Commit message:** `fix(solve): reject RHS callables that return a field on a different grid`

## Group 3 — tests

- New `test/integrator_grid_mismatch_test.cpp`:
  - Throw tests for `integrate(... RK4)` and `integrate(... ExplicitEuler)` on a mismatched-RHS callable. Three flavors of mismatch:
    1. Different total cell count (`psi` is `8×8×8`, RHS returns `4×4×4`).
    2. Same total cell count, different per-axis extents (`psi` is `8×8×8`, RHS returns `4×16×8`).
    3. Same shape, different cell size (`psi.getCellSize()` is `(0.1, 0.1, 0.1)`, RHS returns a field with `(0.2, 0.2, 0.2)`).
  - Throw test for `solve::detail::axpyFresh` directly on two mismatched fields.
  - Per-stage label sweep on RK4: an RHS that returns a correct grid for `k1` but wrong for `k2` (or `k3`, or `k4`) confirms the helper names the offending stage in its message. (Implemented by a stateful lambda that switches behavior on the n-th call.)
  - Smoke test: a correctly-grid-preserving RHS (e.g., `[](const Field<double>& f){ return diff::laplacian(f); }`) integrates without throwing — confirming the precondition does not regress the existing happy path.

**Commit message:** `test(solve): cover RHS grid-mismatch precondition + per-stage labels`

## Group 4 — version bump + CHANGELOG + STATUS

- Bump `CMakeLists.txt` `VERSION` from `0.14.1` to `0.14.2`.
- Update `test/version_test.cpp` to assert `"0.14.2"`.
- Add a `## 0.14.2 — Fix: silent Grid-mismatch on RHS return in time integrators` block to `CHANGELOG.md` (above the `0.14.1` block) with `### Fixed`, `### Tests`, `### Internal` sections.
- Update `specs/STATUS.md`: add a "Phase 14 follow-up: PR #20" entry, bump the **Last updated** line.

**Commit message:** `chore: release 0.14.2 -- silent Grid-mismatch fix in integrators`
