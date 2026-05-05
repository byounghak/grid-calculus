# Validation: Phase 14 hardening — D-grid mismatch + CFL tightening

## Build

- [ ] CI green on Ubuntu GCC.
- [ ] CI green on Ubuntu Clang.
- [ ] CI green on the `clang-debug-nofft` configuration.
- [ ] Warnings clean under the project's strict warning set.
- [ ] `clang-format` and `clang-tidy` pass.
- [ ] `docs-build` job green.
- [ ] `docs-lint` job green. `\since` tags on `cellLaplacian` and the heterogeneous-D `solve::diffuse` overload note the precondition addition.

## Tests

- [ ] `test/fvm_grid_mismatch_test.cpp` exists and covers:
  - [ ] `fvm::cellLaplacian(psi, D)` throws on D-grid mismatch (three flavors: different total size, same total / different per-axis extents, same shape / different cell size).
  - [ ] `fvm::cellLaplacian` throws on `D[i] <= 0` at any cell, with a flat-index message.
  - [ ] `solve::diffuse(psi, D, ...)` throws on D-grid mismatch (three flavors), with a "solve::diffuse" context.
  - [ ] `solve::diffuse(psi, D, ...)` throws on `D[i] <= 0` (now via `cellLaplacian` on the first stage); error message remains sensible.
  - [ ] Happy-path smoke: a correctly-sized D field integrates one step without throwing.
- [ ] CFL tightening coverage:
  - [ ] **High-contrast D**: a `dt` that exceeded `0.5 / (D_max · Σ(1/h²))` under the old bound but is well below `0.5 / (D_max_face · Σ(1/h²))` under the new bound is accepted by `solve::diffuse`.
  - [ ] **Uniform D regression**: bit-identical CFL behavior to the old `D_max` bound (`harmonicMean(D, D) = D`).
  - [ ] **Existing acceptance tests unaffected**: `test/heterogeneous_diffusion_test.cpp` and `test/fvm_diffusion_acceptance_test.cpp` continue to pass with their existing `dt` choices (smooth varying D where `D_max_face ≈ D_max`).
- [ ] Full test suite passes on `clang-debug` (existing 383 tests + new tests, all green).
- [ ] Full test suite passes on `clang-debug-nofft` (existing 240 tests + new tests).

## Documentation

- [ ] `\throws` lines on `fvm::cellLaplacian` and the heterogeneous-D `solve::diffuse` overload list the new preconditions.
- [ ] The "CFL bound" passage on `solve::diffuse`'s Doxygen names the new `D_max_face` and explains that it is strictly tighter than `D_max` for varying `D` and equal for uniform `D`.
- [ ] `CHANGELOG.md` carries a new `0.14.5` block above the existing `0.14.4` block with `### Fixed`, `### Performance`, `### Tests` sections.
- [ ] `specs/STATUS.md` "Last updated" line bumped; "Phase 14 follow-up #5" entry added; PR #22 entry retroactively gains its merge-commit hash.

## Versioning

- [ ] `CMakeLists.txt` project `VERSION` is `0.14.5`.
- [ ] `test/version_test.cpp` asserts `"0.14.5"`.

## Manual verification

- [ ] On a fresh build: `solve::diffuse(psi, D, dt, 1)` with `D` on a different `Grid` than `psi` raises `std::invalid_argument` from `solve::diffuse` (not from inside the integrator) with a context naming "solve::diffuse: D" and the expected vs. actual `(Nx, Ny, Nz, hx, hy, hz)` tuples.
- [ ] On the same build: `fvm::cellLaplacian(psi, D)` called directly (outside `solve::diffuse`) on a mismatched-grid `D` raises with "fvm::cellLaplacian: D" context.
- [ ] On the same build: `fvm::cellLaplacian(psi, D)` with one `D` cell set to `0` raises `std::invalid_argument` naming the offending flat index.
- [ ] On a uniform `D` field, the trajectory of `solve::diffuse` is byte-identical to the trajectory under the previous `D_max`-based bound (regression check; `D_max_face = D_max` for uniform `D`).
- [ ] On a high-contrast `D` field (one cell at 100, neighbors at 1), the maximum stable `dt` under the new bound is approximately `0.5 / (1.98 · Σ(1/h²))` rather than `0.5 / (100 · Σ(1/h²))`.
