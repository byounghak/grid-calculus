# Requirements: Phase 3 — Domain integration (∫ over the grid)

## Roadmap reference

> ## Phase 3 — Domain integration (∫ over the grid)
>
> - **Goal.** A primitive that turns a scalar field into a single number — the integration backbone of functional evaluation.
> - **Deliverables.**
>   - `func/Integrate.hpp` — `integrate(Field<double>) -> double` with deterministic reduction order (pairwise or Kahan).
>   - Test: `∫ 1 dV` returns domain volume; `∫ sin(kx) dV` returns 0 over a period.
> - **Acceptance.** Reduction is deterministic across thread counts (verified by test).

## Decisions made for this feature

- **Reduction strategy: pairwise default + Kahan overload, tag-dispatched.** `integrate(field)` and `integrate(field, Pairwise{})` use pairwise recursive summation (small-block base case ~64 elements). `integrate(field, Kahan{})` uses Neumaier-style compensated summation. The two tags `gridcalc::func::Pairwise` and `gridcalc::func::Kahan` are empty structs declared in `Integrate.hpp`. Pairwise was chosen as default because its rounding-error bound `O(ε log n)` is sufficient for typical functional evaluations and it is cheap (no per-element overhead beyond the recursion). Kahan is exposed as an overload for future callers integrating long sums of widely-varying-magnitude values (anticipated in Phase 11+ functionals); shipping it now closes out the API at Phase 3 and avoids a public-API revision later.
- **API surface: scalar field only at Phase 3.** Signature is `integrate(const Field<double>&, Tag = {}) -> double`. No `Field<Vec3d>` overload (Phase 11+ vector-field functionals add it then), no callable-integrand overload (that is Phase 4's deliverable — `F[ψ] = ∫ f(ψ, ∇ψ, …) dV`).
- **Threading: single-threaded at Phase 3.** OpenMP is not enabled in the build at Phase 3. The acceptance criterion "deterministic across thread counts" is trivially satisfied with one thread. The cross-thread-count determinism *test* ships as a placeholder (`IntegrateTest.DeterministicAcrossInvocations` — verifies bit-equality across two invocations); Phase 20 (performance pass) wires up OpenMP and replaces this placeholder with an environment-varying test. This split was chosen because Phase 3's value is the deterministic-reduction *primitive*, not the threading; pulling Phase 20 work earlier would expand scope without unlocking new functionality.
- **Quadrature rule: periodic midpoint = trapezoidal on equispaced samples.** Discrete formula is `I_h = h_x h_y h_z ∑_{i,j,k} f(i*h_x, j*h_y, k*h_z)`. Exact for samples positioned at the cell origin (no half-cell offset; matches `Grid`'s sampling convention from Phase 1). On periodic analytic integrands the rule is spectrally accurate, which is why the manufactured-solution test below uses an exact-zero target rather than a convergence sweep.
- **Test set: roadmap-mandated + manufactured-solution accuracy + reducer-agreement + determinism placeholder.** Five tests total (see `plan.md` Group 2).
- **Version bump: `0.3.0`** on merge. New public namespace (`gridcalc::func`) and new public symbols (`integrate`, `Pairwise`, `Kahan`); no existing symbol semantics change. SemVer minor bump applies.

## Non-goals

- **No callable-integrand overload `integrate(Field, F)`.** That signature is Phase 4's headline deliverable (`F[ψ] = ∫ f(ψ, ∇ψ, ∇²ψ) dV`) and shipping it now would blur the phase boundary.
- **No `Field<Vec3d>` overload.** Vector-field integration arrives with the higher-rank functionals in Phase 11+.
- **No actual OpenMP parallelism.** Phase 20 owns performance and threading. The determinism contract is asserted by a placeholder test now and exercised properly in Phase 20.
- **No FFT cross-check.** Phase 9 owns spectral verification.
- **No GPU / SIMD intrinsics.** Same — Phase 20.
- **No convergence sweep.** The periodic midpoint rule is spectrally accurate on periodic analytic integrands; a sweep would just verify machine precision at every refinement and produce unstable log-log slopes. The manufactured-solution accuracy test at one fixed `N` is the substantive accuracy check.

## Deferred questions

- **OpenMP enablement and the real cross-thread-count determinism test.** Deferred to Phase 20.
- **`integrate(Field<Vec3d>)` and other rank-extended overloads.** Deferred to whichever Phase 11+ functional first needs them.
- **Non-periodic quadrature rules (Simpson, higher-order Newton–Cotes).** Useful only once Dirichlet/Neumann fields exist (Phase 18); deferred.
- **Reduction-order alternatives (e.g., balanced binary tree with explicit chunk size).** The pairwise recursion's 64-element base case is a reasonable starting point; Phase 20 may revisit if vectorization benchmarks show a different chunk size dominates.
