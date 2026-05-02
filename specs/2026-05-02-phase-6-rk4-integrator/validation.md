# Validation: Phase 6 — RK4 + generic time integrator interface

## Build

- [ ] CI green on Ubuntu GCC and Ubuntu Clang.
- [ ] Warnings clean under `-Wall -Wextra -Wpedantic -Wconversion`.
- [ ] `clang-format` and `clang-tidy` pass (no new regressions vs the codebase-wide hygiene gap noted in Phase 2's STATUS).
- [ ] `gridcalc` target remains INTERFACE — no new source files.
- [ ] Phase 5's `solve::explicitEuler` continues to compile and behave bit-identically (verified by the existing Phase 5 test suite).

## Tests

- [ ] `IntegratorTest.RK4MatchesAnalyticAtFixedN` — at `N=32, dt=0.05, n_steps=100`, RK4 relative max-norm error vs `exp(-3 D T) ψ_0` is `< 1e-2`.
- [ ] `IntegratorTest.CombinedRefinementOrder` — sweep `N ∈ {16, 32, 64}` with `dt = 0.4 * Tag::diffusionCFLLimit / (D * sum(1/h²))`. Both Euler and RK4 show log-log slope in `[1.8, 2.5]` (spatial-dominated combined order).
- [ ] `IntegratorTest.RK4LessErrorThanEulerAtSameDt` — at `N=32` and a shared `dt` within both integrators' CFL bounds, RK4's relative error is at least 4× smaller than Euler's.
- [ ] `IntegratorTest.RK4TimeAccuracyOrder4` — at `N=64` (small `h`, small spatial floor), sweep `dt ∈ {dt0, dt0/2, dt0/4}` with `dt0` chosen so the time-error dominates the spatial floor at the coarsest step. Log-log slope in `[3.5, 4.5]`. **Satisfies the roadmap's `Δt^4` acceptance.**
- [ ] `IntegratorTest.RK4StabilityGap` — `D, dt, h` chosen so `D dt sum(1/h²) ∈ (0.5, 0.6963)`. `solve::diffuse(..., ExplicitEuler{})` throws; `solve::diffuse(..., RK4{})` succeeds and runs to completion.
- [ ] `IntegratorTest.GenericIntegrateDispatchEulerEqualsExplicitEuler` — `solve::integrate(..., ExplicitEuler{})` produces bit-identical output to `solve::explicitEuler(...)`.
- [ ] All Phase 1–5 tests still pass.

## Documentation

- [ ] Every new public symbol (`solve::ExplicitEuler` tag struct, `solve::RK4` tag struct, `solve::integrate` template, the new `diffuse` overload) carries a Doxygen block with `\brief`, `\param`, `\returns`, `\throws`, `\since 0.6.0`. The existing `solve::explicitEuler` retains `\since 0.5.0` and gains a "now a thin wrapper over `integrate(... ExplicitEuler)`" note.
- [ ] User Guide note `docs/user-guide/notes/phase-6-rk4-integrator.md` covering the tag types, the generic `integrate`, the `diffuse` integrator-choice, and a worked example.
- [ ] Developer Note `docs/developer-note/notes/phase-6-rk4-integrator.md` with the five-section structure (Theory · Math derivation · Algorithm · Design decisions · References) per `specs/CLAUDE.md` step 4d. References section non-empty: at least one external citation.
- [ ] `CHANGELOG.md` entry under `0.6.0` lists `solve::ExplicitEuler`, `solve::RK4`, `solve::integrate`, the templated `diffuse`, and the integrator-specific CFL.
- [ ] `STATUS.md` updated: Phase 6 row marked Done; Phase 7 stated as the new Next action.
- [ ] CI `render-docs` job produces fresh aggregate PDFs that include the Phase 6 chapter in both books.

## Performance

- N/A at Phase 6. RK4's per-step allocation cost is documented in the Developer Note; Phase 20 owns the scratch-pool optimization.
