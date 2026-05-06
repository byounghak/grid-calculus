# Plan: Phase 15 follow-up — fix dangling Leaf + perf-fix fused integrate

Four commit-sized groups. Build and tests stay green between every commit.

## Group 1 — Spec files

- `specs/2026-05-05-fix-leaf-dangling-and-fused-integrate-perf/{requirements, plan, validation}.md`.

**Commit message:** `docs: spec files for 0.15.1 audit-driven fixes`

## Group 2 — Reject rvalue `Field<T>` / `Grid` in the ET leaf factories

- `include/gridcalc/tensor/Expressions.hpp`:
  - Add `Leaf(core::Field<T>&&) = delete;` constructor.
  - Add `template <class T> Leaf<T> field(core::Field<T>&&) = delete;` factory.
  - Add `IdentityField(core::Grid&&) = delete;` constructor.
  - Add `IdentityField identityField(core::Grid&&) = delete;` factory.
  - Doxygen note on the surviving lvalue-only `field` and `identityField` factories restating the lifetime contract ("the wrapped field/grid must outlive any expression containing the leaf") and pointing at the deleted rvalue overload.
- New `test/tensor_leaf_lifetime_test.cpp` with:
  - `static_assert(std::is_constructible_v<Leaf<Mat3d>, const Field<Mat3d>&>)` — lvalue still works.
  - `static_assert(!std::is_constructible_v<Leaf<Mat3d>, Field<Mat3d>&&>)` — rvalue rejected at constructor.
  - `static_assert(!std::is_constructible_v<Leaf<Mat3d>, Field<Mat3d>>)` — also rejected for prvalue.
  - Same matrix for `Leaf<double>`, `Leaf<Vec3d>`.
  - Same matrix for `IdentityField` against `const Grid&` (allowed) vs `Grid&&` / `Grid` (rejected).
  - `static_assert(std::is_invocable_v<decltype(&expr::field<Mat3d>), const Field<Mat3d>&>)` — factory accepts lvalue.
  - `static_assert(!std::is_invocable_v<decltype(&expr::field<Mat3d>), Field<Mat3d>&&>)` — factory rejects rvalue.
  - One run-time `TEST` body using the lvalue factories so GoogleTest registers a non-empty fixture (the asserts fire at compile time, but the file needs at least one runtime test for `gtest_discover_tests` to attach a name).
- Add `tensor_leaf_lifetime_test.cpp` to `test/CMakeLists.txt`.

**Commit message:** `fix(tensor): reject rvalue Field/Grid in expr::field / identityField factories`

## Group 3 — Revert streaming integrate; retract perf claim

The original plan was to rewrite the streaming reducers to drop per-cell `% / /` (the auditor's recommendation). That rewrite was implemented locally and benchmarked alongside materialize-then-reduce on `trace(field)` and the elastic-energy integrand at `N ∈ {64, 128}`. The benchmark surfaced that the streaming path is 2-3× **slower** than materialize-then-reduce on Apple ARM with clang -O3 even after the per-cell-division fix, because the materialize path's contiguous `Field<double>` reduces with much better SIMD utilization (see `requirements.md` for the table). Per the auditor's "add a benchmark before shipping the performance claim" recommendation: the benchmark says retract the claim.

- `include/gridcalc/func/Integrate.hpp`:
  - Drop the `detail::pairwiseSumExpr` and `detail::neumaierSumExpr` helpers introduced in 0.15.0.
  - Revert `func::integrate(expr, tag)` to its 0.15.0 form: `return integrate(tensor::expr::materialize(expr), tag);`. The overload still avoids the rank-2 `Field<core::Mat3d>` intermediates (those never enter the picture; the ET nodes evaluate Mat3d arithmetic per cell on the stack); only the rank-0 scalar integrand `Field<double>` is allocated. That is still a meaningful win versus the eager pre-Phase-15 path, just not the "fully fused, no scalar temporary" win the original Doxygen and CHANGELOG language oversold.
  - Update the file-level Doxygen and the per-overload Doxygen on `func::integrate(expr, tag)` to honestly describe what the overload does and links the streaming question to Phase 21 (the dedicated performance pass).
- Existing tests act as the regression gate:
  - `IntegrateExpr_FusedReduction` and `IntegrateExpr_KahanReduction` in `test/tensor_expressions_test.cpp` use `EXPECT_DOUBLE_EQ(via_expr, via_field)`. After the revert the LHS *is* literally `via_field` --- the comparison becomes trivially true. Both keep passing.
  - `ElasticEnergyTest.MatchesViaExpressionLayer` exercises the ET-route integrate against the `func::evaluate` route. Must still pass.

**Commit message:** `perf(func): revert streaming integrate(expr); retract fused-loop perf claim`

## Group 4 — Bookkeeping

- Bump `project(... VERSION 0.15.0)` → `0.15.1` in root `CMakeLists.txt`.
- Update `test/version_test.cpp` to assert `0.15.1`.
- New `## 0.15.1 — Phase 15 audit-driven fixes` block in `CHANGELOG.md`. Calls out both fixes, the audit context, and the local microbenchmark results in the PR description (not the CHANGELOG body, which stays focused on user-visible behavior).
- Refresh `specs/STATUS.md`: bump "Last updated" date, add a new "Phase 15 follow-up #1" paragraph mirroring the Phase 14 follow-up format under "Where the project stands", append PR # to the merged-PR list line. The "Next action" line still points at Phase 16 (this is a follow-up, not a phase advance).

**Commit message:** `chore: bump version to 0.15.1 and refresh STATUS for the audit-driven fixes`

---

After Group 4: `git push -u origin 2026-05-05-fix-leaf-dangling-and-fused-integrate-perf`, open PR with the local microbenchmark numbers in the description, watch CI per `specs/CLAUDE.md` step 6.
