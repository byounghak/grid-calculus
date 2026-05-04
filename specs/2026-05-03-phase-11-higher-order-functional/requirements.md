# Requirements: Phase 11 — Functional with up-to-4th-order gradients

## Roadmap reference

> ## Phase 11 — Functional with up-to-4th-order gradients
>
> - **Goal.** Extend the functional API to densities depending on $\nabla^3$ and $\nabla^4$ terms.
> - **Deliverables.**
>   - Updated `func/Functional.hpp` accepting callables with extended argument lists.
>   - Worked example: phase-field-crystal free energy $F = \int \tfrac{1}{2}\psi[(q_0^2+\nabla^2)^2 - \epsilon]\psi + \tfrac{1}{4}\psi^4 \,dV$.
> - **Acceptance.** PFC functional value matches a hand-computed reference.

## Decisions made for this feature

- **Argument shape: contracted scalars only.** The new derivative arguments enter the callable as `double` scalars — specifically, the biharmonic $\nabla^4\psi$ as a `double` per point, exactly the return type of Phase 8's `diff::biharmonic(psi)`. No `Field<Tensor3>` / `Field<Tensor4>` is introduced; no packed `std::array` of unique multi-index components is exposed. Rank-3 and rank-4 tensor-valued arguments are deferred to Phase 13/14, which is where the rank-6 functional eval (i.e. functional evaluation over vector/tensor data) actually motivates the design choices for a tensor type. This is the smallest API slice that hits the roadmap acceptance bar — PFC's $(q_0^2+\nabla^2)^2$ kernel only needs $\psi$, $\nabla^2\psi$, and $\nabla^4\psi$ (all scalars).
- **No $\nabla^3\psi$ argument is shipped.** $\nabla^3\psi$ has no useful scalar contraction (the natural rank-1 contraction is $\nabla(\nabla^2\psi)$, a `Vec3d`). Rather than add a `Vec3d` 3rd-derivative argument that PFC does not exercise, Phase 11 omits the 3rd-derivative slot entirely. If a future phase needs $\nabla(\nabla^2\psi)$ or a true rank-3 tensor argument, it lands in Phase 13/14 alongside the tensor-type machinery.
- **Arity dispatch: discrete SFINAE arities, extending Phase 4's pattern.** `func::evaluate` gains one new branch — `f(double psi, const Vec3d& grad_psi, double lap_psi, double biharm_psi)` — added at the **top** of the `if constexpr` ladder so the longest-arity branch wins under the existing tie-breaking rule. The four supported arities at Phase 11: `f(ψ)`, `f(ψ, ∇ψ)`, `f(ψ, ∇ψ, ∇²ψ)`, `f(ψ, ∇ψ, ∇²ψ, ∇⁴ψ)`. The variadic / context-object alternative (a single `FieldPoint` struct with lazy accessors) is a larger redesign of Phase 4's surface and was rejected for Phase 11 — the discrete-arity ladder mirrors existing code, materializes only what the matched arity needs, and requires the smallest delta to the static-assert error message.
- **Eager materialization of the biharmonic.** Same model as Phase 4: when the 4-arg branch is matched, `evaluate` calls `diff::biharmonic(psi)` once to materialize a fresh `Field<double>`, then forms the integrand pointwise. Per-point lazy stencils, fused integrand+reducer loops, and contraction expression templates remain Phase 14 deliverables.
- **No contraction-expression-template machinery in Phase 11.** The roadmap's "Decisions worth knowing" entry calls out contraction ETs for "Phases 11 and 14" jointly. Reviewed for Phase 11: with the contracted-scalar-only argument shape, no rank-3 / rank-4 intermediate tensor is ever materialized, so there is no rank-6 explosion to fend off. The ET layer's design choices are driven by the rank-2 (and higher) field arguments that arrive in Phase 13/14 — building it now without those use cases would be premature. No `tensor/Contraction.hpp` lands in this PR.
- **Worked example: phase-field-crystal on `ψ = sin(x) sin(y) sin(z)`** over `[0, 2π]³`, with parameters $q_0 = 1$, $\epsilon = 0.1$. Hand-computed reference:
  - $\int \psi^2\, dV = \pi^3$ (each axis contributes $\int_0^{2\pi}\sin^2 = \pi$).
  - $\nabla^2\psi = -3\psi$ (sum of three $-\sin$ contributions); $\nabla^4\psi = \nabla^2(-3\psi) = 9\psi$.
  - $\int \psi\,\nabla^2\psi\, dV = -3\pi^3$; $\int \psi\,\nabla^4\psi\, dV = 9\pi^3$.
  - $\int \psi^4\, dV = \big(\tfrac{3\pi}{4}\big)^3 = \tfrac{27\pi^3}{64}$, using $\int_0^{2\pi}\sin^4 = \tfrac{3\pi}{4}$.
  - Expanding $(q_0^2+\nabla^2)^2 = q_0^4 + 2q_0^2\nabla^2 + \nabla^4$:
    $$F = \int\!\Big[\tfrac{1}{2}(q_0^4-\epsilon)\psi^2 + q_0^2\,\psi\nabla^2\psi + \tfrac{1}{2}\psi\nabla^4\psi + \tfrac{1}{4}\psi^4\Big]\, dV.$$
  - With $q_0 = 1$, $\epsilon = 0.1$: $F_{\text{ref}} = \big[\tfrac{1}{2}(1 - 0.1) + (-3) + \tfrac{1}{2}\cdot 9 + \tfrac{27}{4 \cdot 64}\big]\,\pi^3 = (0.45 - 3 + 4.5 + \tfrac{27}{256})\,\pi^3 = 2.05546875\,\pi^3 \approx 63.7383$. All four functional terms contribute non-trivially, so the test exercises every term independently.
- **Stencil-order error tolerance.** At `N = 64` with the band-limited ψ, the dominant discrete error in $F$ comes from the biharmonic term (Phase 8's fused 6-term direct stencil at default `Order = 2`, leading-order truncation $O(h^2)$). With $h = 2\pi/64$, $h^2 \approx 9.6\times 10^{-3}$; the leading-error constant on $\nabla^4\sin$ is bounded by $|k|^6 / 12$ (analogous to Phase 8's biharmonic acceptance derivation). This puts the expected relative error in $F$ at the $\sim 10^{-3}$–$10^{-2}$ level. The hand-computed test sets the tolerance at `2e-2` relative; the convergence sweep confirms the $O(h^2)$ rate.
- **Verification depth: hand-computed reference + spectral PFC cross-check (independent numerical path).** Two test patterns:
  1. **`PfcHandComputed`** — `func::evaluate(psi, pfc_lambda)` (production FD path) versus the analytical $F_{\text{ref}}$. Tolerance `2e-2` relative at `N = 64`.
  2. **`PfcSpectralCrossCheck`** — re-materializes the same PFC integrand using `spectral::laplacian(psi)` and `spectral::biharmonic(psi)` (Phase 9), assembles the integrand by hand, integrates with `func::integrate`, and asserts agreement with $F_{\text{ref}}$ at `< 1e-12` relative. Confirms the analytical reference is correct via a second independent numerical path; spectral is exact for band-limited inputs so the residual is at machine-precision integration rounding. The existing FD–FFT cross-check fixture (`test/fd_fft_crosscheck_test.cpp`) is **not** extended — it parameterizes over individual FD operators per the Phase 9 design, and `func::evaluate` is not itself an FD operator. The PFC spectral cross-check lives in the Phase 11 test file.
- **No `Tag` change.** The `Pairwise` / `Kahan` reducer tag and forwarding to `func::integrate` are unchanged from Phase 4; the new arity is appended to the existing `if constexpr` ladder.
- **Doxygen `\since` accounting.** `func::evaluate` itself stays `\since 0.4.0 (function); 0.11.0 (4-arg arity)`. The function-level `\file` tag is updated similarly. No new public symbol is introduced — only the dispatch ladder gains a branch — so there is no fresh `\since 0.11.0` declaration to add. The `\since`-lint check (Phase 10's permanent CI gate) continues to pass without changes to `scripts/check-since.py`.
- **Version bump: `0.11.0`** on merge. New behavior (extended arity); no existing arity semantics change.

## Non-goals

- **No tensor type (`Tensor3` / `Tensor4`).** Deferred to Phase 13/14.
- **No rank-3 / rank-4 callable arguments** (whether as `std::array<double, K>` or as a tensor type). Phase 11 ships only the contracted-scalar 4-arg arity.
- **No $\nabla(\nabla^2\psi)$ (gradient of Laplacian) callable argument.** A natural rank-1 3rd-derivative argument is not motivated by PFC; defer to whichever future phase first needs it.
- **No `tensor/Contraction.hpp`** (single / double contraction operations). Deferred to Phase 14.
- **No expression-template / lazy-evaluation refactor of `func::evaluate`.** Eager materialization remains the model.
- **No vector- or tensor-field functionals.** Deferred to Phase 14.
- **No variational derivative `δF/δψ`.** Deferred to Phase 12.
- **No FD–FFT cross-check fixture extension.** The Phase 9 fixture is operator-level; functional-level cross-checks live in the Phase 11 test file.
- **No CFL or stability work for PFC time evolution.** Phase 11 evaluates the *functional value* only — no time-stepping, no $\delta F/\delta\psi$. PFC dynamics are out of scope until a later worked example.

## Deferred questions

- **Tensor type and tensor-valued callable arguments.** Phase 13 introduces vector/tensor field types; Phase 14 wires them into `func::evaluate`.
- **Contraction expression templates.** Phase 14, alongside the tensor-functional support that motivates them.
- **`∇(∇²ψ)` (rank-1 third-derivative) callable argument.** Deferred to whichever functional first needs it; none of Phases 12–16 are committed to require it.
- **Variadic / context-object callable interface.** Considered for Phase 11 (a single `FieldPoint` struct with lazy `psi() / grad() / lap() / biharm()` accessors) but rejected as a larger redesign of Phase 4's surface than this phase's deliverables warrant. May be reconsidered in Phase 14 once the parameter pack widens further.
- **PFC time evolution and $\delta F/\delta\psi$.** Phase 11 covers the functional value only; dynamic PFC simulations are deferred indefinitely (no committed phase).
