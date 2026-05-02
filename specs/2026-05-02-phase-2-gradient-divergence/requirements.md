# Requirements: Phase 2 — Gradient and divergence (scalar)

## Roadmap reference

> ## Phase 2 — Gradient and divergence (scalar)
>
> - **Goal.** Complete the basic first-order operators.
> - **Deliverables.**
>   - `diff/Gradient.hpp` — `gradient(field) -> Field<Vector3d>`.
>   - `diff/Divergence.hpp` — `divergence(vfield) -> Field<double>`.
>   - Round-trip identity test: `divergence(gradient(ψ)) ≈ laplacian(ψ)` to expected order.
> - **Acceptance.** Convergence-order test passes for both operators on trig inputs.

## Decisions made for this feature

- **Gradient return type is `Field<Vec3d>`.** Single field-of-vectors, matching the roadmap. `Field<T>` is already templated from Phase 1, so this is a one-line return — no wrapper struct, no `getComponent` helper. Three-component decompositions, if ever needed, can land later as additive helpers without changing the canonical API.
- **`Vec3d` lives in `gridcalc/core/EigenAliases.hpp` from Phase 2 forward.** The alias was originally declared in `core/Grid.hpp` (Phase 1) with a comment promising relocation when more operators needed it; Phase 2 is that moment. `core::Vec3d` keeps the same fully-qualified name so Phase 1 source compatibility is preserved. Future vector/tensor typedefs (e.g., `Mat3d`, `Tensor3` if those land in Phase 13/14) belong in this same header.
- **Divergence input is a single `Field<Vec3d>`.** No per-axis-component overload at Phase 2 — symmetric with the gradient return, keeps the API minimal, and is sufficient for the round-trip test. A 3-component overload would just duplicate the implementation; if a downstream caller ever needs it, it can be added then.
- **First-derivative stencil is a new sibling template, `stencil::FirstDerivative<Order>`, in a new header `stencil/FirstDerivative.hpp`.** Phase 1's `stencil::Coefficients<Order>` is documented as the **second-derivative** table; reusing the name for first-derivative coefficients would silently change semantics for existing callers and break Phase 1's public contract as recorded in `STATUS.md`. Two narrowly-named tables alongside each other is the clearer move at this phase. A unifying redesign (e.g., a 2-parameter `stencil::Coefficients<Derivative, Order>`) is deferred to Phase 7, when both tables gain `Order = 4` specializations and the design tradeoffs are easier to evaluate.
- **Convergence-test sweep mirrors Phase 1.** `N ∈ {16, 32, 64}` over a `2π` cube on `ψ = sin(x) sin(y) sin(z)` for the gradient sweep and the round-trip sweep. The independent divergence test uses `V = (sin x cos y, cos x sin y, sin z)` so divergence is exercised on a vector field that is *not* the output of `gradient`.
- **Round-trip identity is verified via convergence to the analytical Laplacian, not the discrete one.** At finite `h`, `D₁ ∘ D₁ ≠ D₂` — the composed operator uses a stride-2 stencil (effective half-width 2, weights `(1, 0, -2, 0, 1)/(4h²)`) while the direct Laplacian uses the 3-point stencil. Both converge to `∇²` at order 2 but with different leading constants. The test therefore compares `divergence(gradient(ψ))` to the analytical `-(kx²+ky²+kz²) ψ`, not to `laplacian(ψ)`. This still answers the roadmap's "to expected order" requirement.
- **Version bump: `0.2.0`** on merge. Phase 1 shipped as `0.1.0`; Phase 2 adds new public headers (`diff/Gradient.hpp`, `diff/Divergence.hpp`, `stencil/FirstDerivative.hpp`, `core/EigenAliases.hpp`) without changing existing semantics, so SemVer minor bump applies.

## Non-goals

- **No `gradient` overload returning three `Field<double>`s.** Callers that need per-axis access can iterate over the `Field<Vec3d>` and project; if a typed component view becomes worth the API surface, add it then.
- **No `divergence` overload for three `Field<double>` arguments.** Same reasoning.
- **No 4th-order accuracy.** `FirstDerivative<4>` is deferred to Phase 7 (which already plans to add `Coefficients<4>`); the primary template stays undefined so unsupported orders trip at compile time at the call site.
- **No FFT/spectral cross-check.** That gate arrives in Phase 9.
- **No benchmarks.** Performance work is centralized in Phase 20.
- **No vector field operator beyond divergence.** Curl and the Hessian land in Phase 5/6 per the roadmap.

## Deferred questions

- **Unified `stencil::Coefficients<Derivative, Order>` redesign.** Sibling tables `Coefficients<Order>` (2nd derivative) and `FirstDerivative<Order>` (1st derivative) are intentionally separate at Phase 2. Phase 7 (4th-order accuracy) is the natural moment to revisit a unifying template.
- **`stencil::FirstDerivative<4>`.** Deferred to Phase 7 alongside `Coefficients<4>`.
- **Per-component access helpers on `Field<Vec3d>` (e.g., `getComponent(field, axis)`).** Deferred until Phase 13 (vector/tensor fields) shows whether functional callers actually need it.
- **`Mat3d` / tensor typedefs in `core/EigenAliases.hpp`.** Added when Phase 13/14 lands, not preemptively.
- **Decoupled-per-sublattice gradient mode.** Phase 16, per `STATUS.md` "Decisions worth knowing."
