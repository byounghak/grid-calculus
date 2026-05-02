# Plan: Phase 6 — RK4 + generic time integrator interface

## Group 1 — Add `solve/Integrator.hpp` (tag types + generic dispatcher)

- Create `include/gridcalc/solve/Integrator.hpp`. Defines two empty tag structs in the `gridcalc::solve` namespace:
  - `struct ExplicitEuler { static constexpr double diffusionCFLLimit = 0.5; };`
  - `struct RK4         { static constexpr double diffusionCFLLimit = 0.6963; };` (the heat-equation stability bound for RK4 — derived in the Developer Note).
- Declare the generic dispatcher
  `template <class Rhs, class Integrator> void integrate(Field<double>& psi, Rhs&& rhs, double dt, int n_steps, Integrator tag);`
  with the actual definitions provided in `solve/ExplicitEuler.hpp` (Group 2 refactor) and `solve/RK4.hpp` (Group 3) as overloads on the tag.
- `Integrator.hpp` `#include`s both `ExplicitEuler.hpp` and `RK4.hpp` so callers only need this one header to get everything.
- Doxygen blocks; `\since 0.6.0`.

**Commit message:** `feat(solve): add Integrator.hpp tag types + generic integrate dispatcher`

## Group 2 — Refactor `solve/ExplicitEuler.hpp`

- Update `include/gridcalc/solve/ExplicitEuler.hpp`:
  - Add overload `template <class Rhs> void integrate(Field<double>& psi, Rhs&& rhs, double dt, int n_steps, ExplicitEuler);` containing the existing per-step loop body. Move the SFINAE check + dt/n_steps validation into this overload.
  - Refactor the existing `explicitEuler(...)` free function to be a thin wrapper: `inline void explicitEuler(...) { integrate(psi, rhs, dt, n_steps, ExplicitEuler{}); }`. Phase 5 callers and the Phase 5 test suite continue to compile unchanged.
- Doxygen on the new `integrate` overload; preserve Phase 5's `\since 0.5.0` on `explicitEuler` and add `\since 0.6.0` on the dispatch overload.
- Header now includes `Integrator.hpp` for the `ExplicitEuler` tag definition.

**Commit message:** `refactor(solve): explicitEuler becomes thin wrapper over integrate(... ExplicitEuler)`

## Group 3 — Add `solve/RK4.hpp`

- Create `include/gridcalc/solve/RK4.hpp`. Defines
  `template <class Rhs> void integrate(Field<double>& psi, Rhs&& rhs, double dt, int n_steps, RK4);`.
- Implementation: classic 4-stage Runge–Kutta:
  ```
  k1 = rhs(psi)
  k2 = rhs(psi + (dt/2) * k1)
  k3 = rhs(psi + (dt/2) * k2)
  k4 = rhs(psi +  dt    * k3)
  psi += (dt/6) * (k1 + 2 k2 + 2 k3 + k4)
  ```
  Each `psi + a*k_i` is materialized as a temporary `Field<double>`; `k_i` is whatever `rhs(...)` returns. Per-step allocation is 4 stages × 1–2 fields ≈ 6 fields. Phase 20 will optimize.
- SFINAE check on the RHS callable signature, mirroring `explicitEuler`'s.
- Doxygen with `\since 0.6.0`.

**Commit message:** `feat(solve): add classic RK4 implementation as integrate(... RK4)`

## Group 4 — Refactor `solve/Diffusion.hpp` to accept an integrator tag

- Update `include/gridcalc/solve/Diffusion.hpp`:
  - Promote `diffuse` to a template:
    `template <class Integrator = ExplicitEuler> void diffuse(Field<double>& psi, double D, double dt, int n_steps, Integrator tag = {});`.
  - The CFL check reads `Integrator::diffusionCFLLimit` rather than the hard-coded `0.5`. The error message names the actual limit and the integrator tag's class name (via a `static_assert` + a `name()` helper, or just by formatting `D*dt*sum(1/h²)` and the per-tag limit).
  - Forward to `integrate(psi, &diff::laplacian, D * dt, n_steps, tag)`.
- Phase 5 callers (no integrator tag) continue working; the default tag is `ExplicitEuler{}`.

**Commit message:** `refactor(solve): diffuse() now templated on integrator tag (default ExplicitEuler)`

## Group 5 — Tests

- Add `test/integrator_test.cpp` (six new tests):
  - `IntegratorTest.RK4MatchesAnalyticAtFixedN` — at `N = 32, D = 0.1, dt = 0.05, n_steps = 100`, RK4 produces relative max-norm error `< 1e-2` against `exp(-3 D T) ψ_0` (much tighter than Phase 5's `5e-2` bound for Euler at the same parameters).
  - `IntegratorTest.CombinedRefinementOrder` — sweep `N ∈ {16, 32, 64}` with `dt = 0.4 * Tag::diffusionCFLLimit / (D * sum(1/h²))` (each integrator at its own CFL-bound fraction). Both Euler and RK4 show log-log slope in `[1.8, 2.5]` (combined error is spatial-dominated at order 2). Verifies the user-requested combined-refinement test.
  - `IntegratorTest.RK4LessErrorThanEulerAtSameDt` — at `N = 32` and a shared `dt` well within Euler's CFL limit, RK4's relative error is at least 4× smaller than Euler's. Demonstrates the time-accuracy benefit even when both integrators converge at the same combined order.
  - `IntegratorTest.RK4TimeAccuracyOrder4` — at `N = 64` (so spatial error is `~h²/12 ≈ 1.6e-3`, dominated by spatial when `dt^4 < 1.6e-3` — but here the test sweeps from large `dt` where time error dominates so we can see the slope before the spatial-error floor kicks in). Sweep `dt ∈ {dt0, dt0/2, dt0/4}` with `dt0` chosen so each level has time-error above the spatial floor; check log-log slope in `[3.5, 4.5]`. **This test satisfies the literal roadmap acceptance "error scales as `Δt^4`."**
  - `IntegratorTest.RK4StabilityGap` — pick `D, dt, h` so `D dt sum(1/h²)` lies between `0.5` (Euler) and `0.6963` (RK4): Euler throws, RK4 runs without throwing. Demonstrates the integrator-specific CFL.
  - `IntegratorTest.GenericIntegrateDispatchEulerEqualsExplicitEuler` — `solve::integrate(psi, rhs, dt, 100, ExplicitEuler{})` produces bit-identical output to `solve::explicitEuler(psi, rhs, dt, 100)` on the same input. Verifies the wrapper is genuinely a no-op.
- Wire `integrator_test.cpp` into `test/CMakeLists.txt`.

**Commit message:** `test(solve): add RK4 + generic integrator dispatch tests`

## Group 6 — Phase 6 doc-notes

- `docs/user-guide/notes/phase-6-rk4-integrator.md` — user-facing summary: tag types, signature of `solve::integrate`, the choice between Euler and RK4 in `solve::diffuse`, a worked example using RK4 to integrate the heat equation. Includes a brief explanation of why RK4 is preferable when accuracy matters (4× cost per step buys 4 orders of accuracy reduction).
- `docs/developer-note/notes/phase-6-rk4-integrator.md` — five-section structure per `specs/CLAUDE.md` step 4d:
  1. **Theory**: Runge–Kutta family; classic RK4 Butcher tableau; why intermediate-stage evaluations of `rhs` improve the temporal Taylor match.
  2. **Math derivation**: explicit Taylor expansion of the RK4 update, showing match to $\dot\psi$ through $\Delta t^4$; derivation of RK4's stability region intersection with the negative real axis at `~-2.7853` (and hence the heat-equation CFL bound `D dt sum(1/h_a^2) <= 2.7853 / 4 ≈ 0.6963`); discussion of why combined refinement at CFL-limited `dt` is spatial-dominated for both integrators.
  3. **Algorithm**: the four-stage update, intermediate-field allocation cost, the tag-dispatched `solve::integrate` overload structure, the `Tag::diffusionCFLLimit` constexpr.
  4. **Design decisions**: tag dispatch over class hierarchy (cite `requirements.md`); `explicitEuler` kept as thin wrapper; `diffuse` templated rather than two drivers; combined-refinement convergence test plus a separate dt-only sweep to satisfy the roadmap's literal `Δt^4` acceptance.
  5. **References**: Hairer, Nørsett & Wanner *Solving Ordinary Differential Equations I* (textbook citation, Chapter II.1 on classical RK methods); Butcher *Numerical Methods for Ordinary Differential Equations* (textbook citation); a permanent-URL reference for RK4's stability region (Wikipedia or a Wolfram MathWorld page).

**Commit message:** `docs(notes): add Phase 6 User Guide and Developer Note notes`

## Group 7 — Bookkeeping

- Bump `CMakeLists.txt` `project(gridcalc VERSION 0.5.0)` → `0.6.0`.
- Update `test/version_test.cpp` (`ReturnsZeroFiveZero` → `ReturnsZeroSixZero`).
- `CHANGELOG.md`: new `## 0.6.0 — Phase 6: RK4 + generic time integrator interface` entry.
- `specs/STATUS.md`: bump "Last updated"; move Phase 6 row to Done; update "Where the project stands" + "Next action" (Phase 7 — higher-order accuracy stencils).

**Commit message:** `chore: bump version to 0.6.0 and refresh STATUS for Phase 6`
