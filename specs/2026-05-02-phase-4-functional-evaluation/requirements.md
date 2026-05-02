# Requirements: Phase 4 — Functional evaluation, callable API (scalar, periodic)

## Roadmap reference

> ## Phase 4 — Functional evaluation, callable API (scalar, periodic)
>
> - **Goal.** First end-to-end functional: $F[\psi] = \int f(\psi, \nabla\psi, \nabla^2\psi)\,dV$ with user-supplied $f$.
> - **Deliverables.**
>   - `func/Functional.hpp` — accepts a callable `f(double, Vector3d, double)` and evaluates pointwise then integrates.
>   - Worked example: Ginzburg–Landau free energy $F = \int \tfrac{1}{2}|\nabla\psi|^2 + W(\psi^2-1)^2 \,dV$.
>   - Tutorial doc page walking through the example.
> - **Acceptance.** Numerical $F$ matches a hand-computed reference for a chosen $\psi$ field within stencil-order error.

## Decisions made for this feature

- **Function name: `gridcalc::func::evaluate`.** Verb-style, matches the project's `getX`/`isY`/`buildZ` convention (CLAUDE.md ref). Pairs naturally with `gridcalc::func::integrate` from Phase 3 — `evaluate` is the higher-level operation whose final step is calling `integrate`.
- **SFINAE-detected callable arity.** `evaluate` accepts callables with three signatures: `f(double)`, `f(double, const Vec3d&)`, `f(double, const Vec3d&, double)`. The implementation uses `if constexpr` + `std::is_invocable_r_v<double, F&, …>` to pick the matching dispatch. This is a strict superset of the roadmap's "callable `f(double, Vector3d, double)`" — it gives users the semantically natural signature for low-order functionals (e.g., `f(ψ)` for `∫ ψ² dV`) and skips the corresponding stencil materialization. Invocations that don't match any of the three signatures trip a clean compile-time error via a `detail::deferred_false<F>` `static_assert` (the deferred-false trick is required because plain `static_assert(false, …)` inside `if constexpr` fires unconditionally under C++17).
- **Vec3d argument is `const Vec3d&`.** The user's lambda may declare it as `const Vec3d&` or by-value `Vec3d` — both work because `const Vec3d&` binds to either. The SFINAE check is performed against `const Vec3d&` so the arity dispatch is unambiguous.
- **Eager materialization.** When `evaluate` needs `∇ψ`, it internally calls `diff::gradient(psi)` to materialize a `Field<Vec3d>`; when it needs `∇²ψ`, it internally calls `diff::laplacian(psi)` to materialize a `Field<double>`. The integrand `f(...)` is then computed pointwise into a fresh `Field<double>` and passed to `func::integrate`. This reuses Phase 1–3 primitives directly and keeps the implementation under ~50 lines. Per-point lazy stencil application would save the two intermediate field allocations but at the cost of duplicating the gradient/Laplacian formulas inline; the optimization is deferred to Phase 11+ when expression-template plumbing exists.
- **Reduction tag is forwarded to `integrate`.** `evaluate` accepts `Tag tag = Pairwise{}` (the existing tag types from `func/Integrate.hpp`) and forwards it to the final `integrate` call. A `static_assert` constrains `Tag` to `Pairwise` or `Kahan` so nonsense tag types fail at compile time at the call site.
- **Worked example: Ginzburg–Landau on `ψ = sin(x) sin(y) sin(z)`.** The roadmap deliverable. Hand-computed reference (with `W = 1`): $F = \tfrac{3}{2}\pi^3 + \tfrac{411}{64}\pi^3 = \tfrac{507}{64}\pi^3 \approx 245.62$. Derivation:
  - $\int (\tfrac{1}{2})|\nabla\psi|^2\, dV = \tfrac{1}{2} \cdot 3\pi^3 = \tfrac{3}{2}\pi^3$ (each axis-aligned $\cos^2 \cdot \sin^2 \cdot \sin^2$ term contributes $\pi^3$, three axes total).
  - $\int W(\psi^2-1)^2\, dV = W\cdot[\int \psi^4 - 2\int\psi^2 + \int 1] = W\cdot[\tfrac{27}{64}\pi^3 - 2\pi^3 + 8\pi^3] = W\cdot\tfrac{411}{64}\pi^3$, using $\int_0^{2\pi}\sin^4 = \tfrac{3\pi}{4}$, $\int_0^{2\pi}\sin^2 = \pi$, $\int_0^{2\pi}1 = 2\pi$.
- **Test set: GL hand-computed at one fixed N + 2nd-order convergence sweep + per-arity sanity tests + reducer agreement.** Six tests total (see `plan.md` Group 2).
- **Stencil-order error tolerance.** At `N = 32` with `ψ = sin(x) sin(y) sin(z)`, the dominant discrete error in the GL functional is from the squared-gradient term: leading-order relative error in `|∇ψ|²` is `~h²/3` per axis (the gradient component error `-k³h²/6` enters squared), so the integrated `(1/2)|∇ψ|²` term has relative error `~h²/3 ≈ 1.3%` at `N=32`. The potential term integrates spectrally. The hand-computed test sets the tolerance at `1e-2` relative, which is met with headroom; the convergence sweep confirms the `O(h²)` rate.
- **Version bump: `0.4.0`** on merge. New public symbol (`func::evaluate`); no existing symbol semantics change.

## Non-goals

- **No expression-template / lazy-evaluation optimization.** Eager materialization is the Phase 4 design; per-point lazy stencils, fused integrand+reducer loops, and contraction expression templates are all Phase 11+ deliverables.
- **No vector-field functionals (`F[V] = ∫ f(V, ∇V, …) dV`).** Phase 11+ extends the API to vector and tensor fields once those operators land.
- **No higher-order derivatives in `f`.** Phase 4 supports `(ψ, ∇ψ, ∇²ψ)` only. Phase 11 adds up-to-4th-order gradients; Phase 6 adds biharmonic and mixed partials.
- **No automatic differentiation / variational derivative `δF/δψ`.** Functional derivatives are out of scope until Phase 12 (CH demo) needs them, where they will be implemented from explicit closed-form expressions, not by AD.
- **No threading.** Phase 20 owns performance; the integrator's reduction inherits the Phase 3 single-threaded behavior and the same determinism placeholder.

## Deferred questions

- **Per-point lazy / fused-loop optimization.** Deferred to Phase 11+ where expression-template machinery is being built anyway.
- **`evaluate(Field<Vec3d>, …)` for vector fields.** Deferred to whichever Phase 11+ functional first needs it.
- **Variational derivatives `δF/δψ`.** Deferred to Phase 12 (CH demo).
- **Higher-order callable arities (e.g., `f(ψ, ∇ψ, ∇²ψ, ∇³ψ)`).** Deferred to Phase 11.
