# Validation: Phase 15 follow-up — fix dangling Leaf + perf-fix fused integrate

## Build

- [ ] CI green on Ubuntu GCC and Ubuntu Clang.
- [ ] `docs-build` and `docs-lint` jobs green; the Doxygen note on the surviving lvalue-only `field` / `identityField` factories renders cleanly.
- [ ] `\since`-tag lint passes — no new public symbols (the `= delete`d overloads are surface restrictions, not new symbols); existing symbols keep their `\since 0.15.0`.
- [ ] Phase 9's FD–FFT cross-check fixture still passes — this PR adds no new FD operators.
- [ ] `GRIDCALC_BUILD_EXAMPLES=ON` builds `elastic_energy` cleanly under the project's strict warning set on both CI compilers.
- [ ] Warnings clean under `-Wall -Wextra -Wpedantic -Wconversion`.
- [ ] `clang-format` / `clang-tidy` clean.

## Tests

### Lifetime contract (Group 2)

- [ ] `static_assert` matrix in `test/tensor_leaf_lifetime_test.cpp` covering:
  - `Leaf<T>` for `T ∈ {double, Vec3d, Mat3d}`: `const Field<T>&` constructible (true), `Field<T>&&` and prvalue `Field<T>` not constructible (false).
  - `IdentityField`: `const Grid&` constructible (true), `Grid&&` and prvalue `Grid` not constructible (false).
  - `expr::field<T>` factory: invocable on `const Field<T>&` (true), not invocable on `Field<T>&&` (false).
  - `expr::identityField` factory: invocable on `const Grid&` (true), not invocable on `Grid&&` (false).
- [ ] One run-time `TEST` body uses the lvalue factories so the file shows up in `gtest_discover_tests`.

### Reverted streaming integrate (Group 3)

- [ ] `IntegrateExpr_FusedReduction` (in `test/tensor_expressions_test.cpp`, pre-existing) — `func::integrate(0.5 * doubleContract(leafA, leafB))` matches `func::integrate(materialize(...))` to `EXPECT_DOUBLE_EQ`. After the revert the LHS *is* `func::integrate(materialize(...))`; the equality is trivially true and the test continues to pass.
- [ ] `IntegrateExpr_KahanReduction` (pre-existing) — same with `Kahan` tag.
- [ ] `ElasticEnergyTest.MatchesViaExpressionLayer` (pre-existing) — `func::evaluate(u, energy_density)` route matches the ET-`func::integrate(0.5 lambda trace²(eps) + mu eps:eps)` route to round-off. The ET route now goes through `materialize` internally; the round-trip still holds.
- [ ] All other Phase 15 ET tests remain green (17 in `tensor_expressions_test.cpp` + 8 in `elastic_energy_test.cpp` + 1 in `elastic_energy_demo_test.cpp` + 6 in `func_evaluate_vector_test.cpp` + 17 in `tensor_grid_mismatch_test.cpp`).

### Regression

- [ ] All previously-passing tests remain green (446 prior on `clang-debug`, 303 prior on `clang-debug-nofft`). This PR adds 1 runtime test (the lifetime fixture) and 16 compile-time `static_assert`s; it removes 2 internal helpers (`detail::pairwiseSumExpr`, `detail::neumaierSumExpr`) that were never user-facing.

## Performance (PR description, this `requirements.md`, and CHANGELOG)

- [ ] Local microbenchmark table recorded in `requirements.md` comparing streaming (post-rewrite) vs.\ materialize-then-integrate on:
  - Cheap integrand: `func::integrate(trace(field))` on a `Field<Mat3d>` at `N ∈ {64, 128}`.
  - Elastic-energy integrand: `func::integrate(0.5 lambda * trace(eps) * trace(eps) + mu * doubleContract(eps, eps))` at `N ∈ {64, 128}`.
- [ ] CHANGELOG `## 0.15.1` block records the perf-claim retraction with the table headline ("streaming 2-3× slower than materialize on Apple ARM clang -O3").
- [ ] `STATUS.md` "Open / deferred items" gains an entry pointing the streaming-with-vectorization question to Phase 21.
- [ ] No permanent benchmark target committed — Phase 21 owns that surface.

## Documentation

- [ ] Doxygen note on the surviving lvalue-only `field` and `identityField` factories restating the lifetime contract.
- [ ] CHANGELOG `## 0.15.1` block added; calls out both fixes and the audit context.
- [ ] `STATUS.md` updated: "Last updated" bumped to `2026-05-05`; new "Phase 15 follow-up #1" paragraph added under "Where the project stands"; PR list extended.

## Bookkeeping

- [ ] `project(... VERSION 0.15.1)` in root `CMakeLists.txt`.
- [ ] `test/version_test.cpp` asserts `0.15.1`.
- [ ] PR opened, every CI check green.
