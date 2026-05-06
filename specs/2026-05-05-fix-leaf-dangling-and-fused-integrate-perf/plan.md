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

## Group 3 — Nested-loop fused reducers

- `include/gridcalc/func/Integrate.hpp`:
  - `detail::pairwiseSumExpr` base case: compute `(i_start, j_start, k_start)` from `start` *once* via the existing `% / /` arithmetic, then walk the next `n` cells incrementing `i` with rollover into `j` and `k`. Recursion structure (split at `n / 2`, base-case threshold `kBaseCase = 64`) unchanged.
  - `detail::neumaierSumExpr` becomes a triple-nested loop `for (int k = 0; k < Nz; ++k) for (int j = 0; j < Ny; ++j) for (int i = 0; i < Nx; ++i)`, evaluating `expr.evalAt(i, j, k)` per cell. Loop signature drops the `(n, Nx, Ny)` arguments in favor of `(grid)` (or just `(Nx, Ny, Nz)`), since we now iterate naturally rather than via flat-index arithmetic.
  - The `func::integrate(expr, tag)` overload calls the rewritten reducers; the tag-dispatch ladder is unchanged.
  - Doxygen on both reducers updated to call out the no-per-cell-division loop pattern and the bit-identical-output contract vs.\ the linear-walk version.
- Existing tests act as the regression gate:
  - `IntegrateExpr_FusedReduction` and `IntegrateExpr_KahanReduction` in `test/tensor_expressions_test.cpp` use `EXPECT_DOUBLE_EQ(via_expr, via_field)`. The rewrite must keep both passing — same order of additions per the i-fastest layout.
  - `ElasticEnergyTest.MatchesViaExpressionLayer` exercises the fused integrate on the elastic-energy integrand and compares to the `func::evaluate` route. Must still pass.

**Commit message:** `perf(func): nested-loop fused integrate without per-cell integer division`

## Group 4 — Bookkeeping

- Bump `project(... VERSION 0.15.0)` → `0.15.1` in root `CMakeLists.txt`.
- Update `test/version_test.cpp` to assert `0.15.1`.
- New `## 0.15.1 — Phase 15 audit-driven fixes` block in `CHANGELOG.md`. Calls out both fixes, the audit context, and the local microbenchmark results in the PR description (not the CHANGELOG body, which stays focused on user-visible behavior).
- Refresh `specs/STATUS.md`: bump "Last updated" date, add a new "Phase 15 follow-up #1" paragraph mirroring the Phase 14 follow-up format under "Where the project stands", append PR # to the merged-PR list line. The "Next action" line still points at Phase 16 (this is a follow-up, not a phase advance).

**Commit message:** `chore: bump version to 0.15.1 and refresh STATUS for the audit-driven fixes`

---

After Group 4: `git push -u origin 2026-05-05-fix-leaf-dangling-and-fused-integrate-perf`, open PR with the local microbenchmark numbers in the description, watch CI per `specs/CLAUDE.md` step 6.
