# Requirements: Phase 6 — RK4 + generic time integrator interface

## Roadmap reference

> ## Phase 6 — RK4 + generic time integrator interface
>
> - **Goal.** Cleaner API for time-stepping; higher-order time accuracy.
> - **Deliverables.**
>   - `solve/Integrator.hpp` interface and `RK4` implementation.
>   - Refactor `Diffusion` to consume the integrator interface.
> - **Acceptance.** Convergence-in-time test: error scales as $\Delta t^4$.

## Decisions made for this feature

- **Tag-dispatched free function.** `solve::integrate(psi, rhs, dt, n_steps, IntegratorTag{})` with empty tag structs `solve::ExplicitEuler` and `solve::RK4`. Same pattern as Phase 3's `Pairwise`/`Kahan` for `func::integrate` and Phase 4's `evaluate`. Compile-time dispatch; no virtual-call overhead; the tag carries integrator-specific properties (e.g., `static constexpr double diffusionCFLLimit`).
- **Phase 5's `solve::explicitEuler` is preserved as a thin wrapper.** Defined in `solve/ExplicitEuler.hpp` as `inline void explicitEuler(psi, rhs, dt, n_steps) { integrate(psi, rhs, dt, n_steps, ExplicitEuler{}); }`. Phase 5 callers and the Phase 5 test suite continue to compile unchanged. Both forms are documented as equivalent — pick whichever reads better at the call site.
- **`solve::diffuse` is templated on the integrator tag, defaulted to `ExplicitEuler`.** Signature becomes `template <class Integrator = ExplicitEuler> void diffuse(Field<double>& psi, double D, double dt, int n_steps, Integrator tag = {})`. Phase 5 calls without a tag retain Phase 5 behavior. The CFL check reads `Integrator::diffusionCFLLimit` so RK4 callers get the larger stability bound automatically.
- **Convergence test setup: combined refinement of `N` and `dt` plus a separate dt-only sweep.** The user-chosen test (combined refinement keeping `T` fixed; CFL-limited `dt` per integrator) demonstrates the **combined order** of accuracy — both Euler and RK4 show slope `~2` because the spatial error dominates at CFL-limited `dt`. To honor the roadmap's literal acceptance bar (`error scales as Δt^4`), a separate **dt-only** sweep at fixed `N` is included that varies `dt` over a window where time error dominates over the spatial-error floor; this exhibits RK4's slope-4. Both tests are part of the validation checklist.
- **RK4 implementation: classic 4-stage Runge–Kutta.** Butcher tableau weights `(1/6, 1/3, 1/3, 1/6)`, stage offsets `(0, dt/2, dt/2, dt)`. Per-step cost: 4 RHS evaluations + 4 intermediate `Field<double>` materializations + a final weighted-sum write to `psi`. Total intermediate allocation per step: `~6` fields (4 stage outputs + 2 intermediate states). Phase 20 owns the optimization to a fixed scratch pool.
- **RK4 heat-equation CFL: `D · dt · sum_a(1/h_a^2) ≤ 0.6963`.** Derived from the RK4 stability region intersecting the negative real axis at `≈ −2.7853`; combined with the discrete Laplacian's `max|λ_h| = 4·sum_a(1/h_a^2)` gives `4 · D·dt · sum(1/h²) ≤ 2.7853`, i.e., the per-tag limit. About `1.39×` looser than Euler's `0.5`. Documented in the Developer Note.
- **Backward compatibility.** Phase 5 tests (`DiffusionTest.GenericExplicitEulerOnZeroRhs`, etc.) must continue passing without modification — they exercise `solve::explicitEuler` directly, which is the wrapper.
- **Version bump: `0.6.0`** on merge. New public symbols (`solve::ExplicitEuler` tag struct, `solve::RK4` tag struct, `solve::integrate` template, `solve::diffuse` newly templated). The `solve::ExplicitEuler` tag *struct* coexists with the existing `solve::explicitEuler` *free function* — same namespace but different identifiers (`ExplicitEuler` vs `explicitEuler`); no clash.

## Non-goals

- **No adaptive time stepping.** Embedded-Pair RK methods (RK4(5), Dormand–Prince) are out of scope. Add only if a downstream phase actually needs error-controlled stepping.
- **No symplectic / structure-preserving integrators.** Heat equation is dissipative; symplectic schemes are for Hamiltonian systems and arrive only if a phase-field flow phase needs them.
- **No multi-step methods (Adams–Bashforth, BDF).** Phase 19 introduces implicit BDF for stiff problems; multi-step explicit methods are not on the roadmap.
- **No virtual / runtime-dispatched integrator interface.** Compile-time dispatch is sufficient; runtime dispatch would force allocation + virtual calls per step.
- **No threading / fused-loop optimization.** Phase 20.

## Deferred questions

- **Persistent scratch buffer for RK4 stages.** Currently each step allocates ~6 fields (one per stage output + 2 for intermediate `psi + a*k_i` evaluations). Phase 20 will reuse a single scratch pool.
- **Exception-safe semantics if `rhs` throws mid-step.** RK4 does 4 RHS evaluations per step; if the third throws, `psi` may be in an indeterminate state (currently it's the value at the start of the step because we don't write back until the final accumulation). Document but don't add explicit rollback.
- **Embedded RK4(5) for adaptive stepping.** Defer until adaptive stepping is actually needed.
- **Interaction with Phase 18 boundary conditions.** Same deferred-test note as Phase 5: RK4 with non-periodic boundary conditions needs its own test set when Phase 18 lands.
