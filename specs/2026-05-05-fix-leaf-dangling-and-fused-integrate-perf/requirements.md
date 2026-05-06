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

- **Streaming `func::integrate(expr, tag)` is reverted to its 0.15.0 form** (forward to `func::integrate(materialize(expr), tag)`). The `detail::pairwiseSumExpr` and `detail::neumaierSumExpr` helpers introduced in 0.15.0 are removed. Reasoning is benchmark-driven: a microbenchmark prototyped per the auditor's recommendation found the streaming path 2-3x **slower** than materialize-then-reduce on the elastic-energy integrand at `N = 64` and `N = 128`, on Apple ARM with clang -O3. Numbers (10-iteration averages, `-O3 -DNDEBUG`):

  | integrand      | grid | streaming (post-rewrite) | materialize+pairwise | streaming/materialize ratio |
  | -------------- | ---- | ------------------------ | -------------------- | --------------------------- |
  | `trace(field)` | 64³  | 1.60 ms                  | 0.74 ms              | 2.2× slower                 |
  | `trace(field)` | 128³ | 6.32 ms                  | 3.67 ms              | 1.7× slower                 |
  | elastic energy | 64³  | 3.69 ms                  | 1.24 ms              | 3.0× slower                 |
  | elastic energy | 128³ | 25.59 ms                 | 10.14 ms             | 2.5× slower                 |

  The materialize path benefits from SIMD-vectorized pairwise reduction over a contiguous `double` array; the streaming path's evalAt-interleaved-with-recursion structure does not vectorize as well under clang -O3. The right fix is *not* to micro-optimize the streaming path further but to reconsider whether streaming is the right approach at all on common hardware. Per the auditor's "add a benchmark before shipping the performance claim" recommendation: the benchmark says retract the claim.

- **Streaming-with-vectorization deferred to Phase 21.** A future re-evaluation may find that an OpenMP-parallelized or hand-SIMD streaming path can beat the contiguous-array materialize path, but exploring that belongs in Phase 21 (the dedicated performance pass), not in a 0.15.1 audit follow-up. The "Open / deferred items" list in `STATUS.md` gains a new entry pinning the question.

- **0.15.0 fused-loop perf claim is retracted in the Doxygen + CHANGELOG.** The Phase 15 chapter and 0.15.0 CHANGELOG block called the ET-integrate overload "fused" and described it as streaming `evalAt` straight into the reducer with "no scalar temporary." That language was at best aspirational and at worst misleading once benchmarks landed. The 0.15.1 CHANGELOG block calls out the retraction; the Doxygen block on `func::integrate(expr, tag)` is rewritten to honestly describe what it does (avoids rank-2 intermediates; allocates one scalar `Field<double>`) and to link the streaming question to Phase 21.

- **No permanent benchmark suite added in this PR.** Phase 21's performance pass owns the benchmark infrastructure (`bench/baselines/` per `tech-stack.md` validation rule). The local microbenchmark above is a throwaway file under `/tmp/bench_integrate_expr.cpp`; only the resulting numbers and decision land in this PR (in this `requirements.md`, the CHANGELOG, and the Doxygen note).

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
