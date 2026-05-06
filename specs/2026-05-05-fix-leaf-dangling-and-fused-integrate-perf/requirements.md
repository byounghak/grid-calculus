# Requirements: Phase 15 follow-up — fix dangling Leaf + perf-fix fused integrate

## Audit reference

Two findings surfaced by an external code audit on the merged Phase 15 surface (PR #24, version `0.15.0`, head commit `d129fbf`):

1. **[high] Lazy field leaves can dangle when wrapping temporaries.** `tensor::expr::Leaf<T>` stores a raw `const Field<T>*` and `expr::field(F)` accepts any `const Field<T>&`, including bound rvalues. Patterns the User Guide chapter explicitly encourages — e.g. `expr::sym(expr::field(diff::gradient(u)))` — wrap the temporary `Field<Mat3d>` returned by `diff::gradient`, destroy it at the end of the full-expression statement, and leave the expression tree pointing into freed storage. Subsequent `materialize` / `func::integrate` walks dereference the dangling pointer (silently wrong values, or crash). The same hazard applies to `IdentityField`, which holds a raw `const Grid*` and accepts an rvalue `Grid` via `expr::identityField(...)`.

2. **[medium] Fused integration pays multiple integer divisions per cell.** Both `detail::pairwiseSumExpr` (base case) and `detail::neumaierSumExpr` reconstruct `(i, j, k)` from a flat linear index using two `%` and two `/` operations *per evaluated cell*. For cheap expressions (e.g. `trace(field)`) this division overhead can dominate the actual reduction work, undermining the "fused-loop" performance claim documented in the Phase 15 chapter and CHANGELOG.

## Decisions made for this feature

- **Reject rvalue `Field<T>` and `Grid` at both the constructor and the factory.**
  - `tensor::expr::Leaf<T>(core::Field<T>&&) = delete;` — constructor.
  - `template <class T> Leaf<T> tensor::expr::field(core::Field<T>&&) = delete;` — factory.
  - `tensor::expr::IdentityField(core::Grid&&) = delete;` — constructor.
  - `IdentityField tensor::expr::identityField(core::Grid&&) = delete;` — factory.

  Defense in depth: deleting at the factory alone closes the documented call site; deleting at the constructor closes the path through any sufficiently determined caller. Both are necessary because `Leaf<T>{rvalue_field}` would otherwise still compile (the constructor is `explicit` but accepts a `const Field<T>&` which an rvalue binds to).

- **Compile-time-only verification.** Runtime tests are impossible by design — the failure is a `delete`d overload, surfacing at compile time. The new test file `test/tensor_leaf_lifetime_test.cpp` uses `static_assert(!std::is_constructible_v<...>)` and `static_assert(!std::is_invocable_v<...>)` to lock the behavior. A single `static_assert(std::is_constructible_v<Leaf<Mat3d>, const Field<Mat3d>&>)` confirms the lvalue path still works.

- **Doxygen restating the lifetime contract.** The surviving lvalue-only factories (`expr::field` and `expr::identityField`) gain an explicit "the wrapped field/grid must outlive the expression" line in their Doxygen blocks, plus a sibling-comment noting that rvalue overloads are deleted.

- **Nested-loop fused reducers, bit-identical to the linear-walk semantics.**
  - `detail::pairwiseSumExpr` retains its linear-range recursion (so the pairwise tree shape — `kBaseCase = 64`, split at `n / 2` — is unchanged). At the base case, compute the starting `(i_start, j_start, k_start)` *once* from `start` via the existing `% / /` arithmetic, then walk the next `n` cells by incrementing `i` with rollover into `j` and `k`. Cost: one pair of `% / /` per ≤64 cells (was: per cell).
  - `detail::neumaierSumExpr` becomes a triple-nested `for (k) for (j) for (i)` loop in i-fastest order. Bit-identical to the current linear walk because `Field<T>::data()` storage and the `linear → (i, j, k)` mapping both use i-fastest order.

  Both rewrites preserve the existing `IntegrateExpr_FusedReduction` and `IntegrateExpr_KahanReduction` `EXPECT_DOUBLE_EQ` regression tests (the order of additions is bit-for-bit identical to the materialize-then-pairwise/Neumaier path).

- **No permanent benchmark suite added in this PR.** Phase 21's performance pass owns the benchmark infrastructure (`bench/baselines/` per `tech-stack.md` validation rule). The PR description records local microbenchmark numbers comparing fused vs.\ materialize-then-integrate on `trace(field)` and the elastic-energy integrand at `N ∈ {64, 128}`, but no `bench/` target lands here. Adding one now would either be tossed or re-implemented when Phase 21 stands up the proper suite.

## Non-goals

- **No new ET node types.** Only the existing `Leaf` / `IdentityField` and the two `detail::*SumExpr` reducers change.
- **No owning-leaf alternative.** The reviewer suggested an "owning leaf" as an alternative if temporary composition is required. Rejected: temporary composition is already explicitly disallowed by the lifetime contract; an owning variant would silently double the storage of every materialized field that flows through the ET layer, which is an antipattern for the "no rank-2 intermediates" guarantee Phase 15's CHANGELOG advertises.
- **No public API change to `func::integrate(expr, tag)`.** Same overload signatures, same numerical output to round-off. The fix is entirely internal to `detail::pairwiseSumExpr` and `detail::neumaierSumExpr`.
- **No widening to other types.** The dangling-Leaf hazard is structural to "raw pointer to a Field"; no other expression node holds raw pointers (combinator nodes own their operands by value), so no other node needs the rvalue rejection.

## Deferred questions

- **Permanent benchmark target for `func::integrate(expr)`.** Pinned to Phase 21.
- **Owning-leaf alternative for temporary composition.** If a future use case actually requires composing temporaries, the right add is a separate `expr::owned(Field<T>)` factory that takes by value and stores by value; not motivated by any current code path.
- **`expr::field` / `expr::identityField` on r-value `Field<T>` / `Grid` const-prvalues from `consteval` / `constexpr` contexts.** Not exercised today.

## Branch and version

- **Branch:** `2026-05-05-fix-leaf-dangling-and-fused-integrate-perf`.
- **Version bump:** `0.15.0 → 0.15.1`.
- **CHANGELOG entry:** new `## 0.15.1 — Phase 15 audit-driven fixes` block.
- **`\since` tag:** the deleted overloads do not need `\since` tags (they are surface restrictions, not new public symbols). Existing public symbols keep their `\since 0.15.0`.
