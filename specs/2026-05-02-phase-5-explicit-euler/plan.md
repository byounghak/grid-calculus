# Plan: Phase 5 — Explicit Euler diffusion solver

## Group 1 — Add `solve::explicitEuler` (generic integrator)

- Create `include/gridcalc/solve/ExplicitEuler.hpp` (new `gridcalc::solve` namespace; new `solve/` subdirectory under the existing public-include tree).
- Define `template <class Rhs> void explicitEuler(core::Field<double>& psi, Rhs&& rhs, double dt, int n_steps)`.
  - Validates `n_steps >= 0` and `dt >= 0`; throws `std::invalid_argument` on violation with a message that includes the offending value.
  - SFINAE check via `std::is_invocable_r_v<core::Field<double>, Rhs&, const core::Field<double>&>`; on failure, `static_assert` via the `detail::deferred_false<F>` trick (same pattern as Phase 4 `evaluate`).
  - Each step: compute `tmp = rhs(psi)`, then `psi[n] += dt * tmp[n]` element-wise via raw pointer access. Per-step allocation cost is one fresh `Field<double>` (Phase 20 will optimize).
- Doxygen blocks on every public symbol with `\since 0.5.0`.

**Commit message:** `feat(solve): add generic explicitEuler integrator`

## Group 2 — Add `solve::diffuse` (diffusion driver) + CFL check

- Create `include/gridcalc/solve/Diffusion.hpp`.
- Define `inline void diffuse(core::Field<double>& psi, double D, double dt, int n_steps)`:
  - Validates `D >= 0`.
  - Calls `detail::checkExplicitDiffusionCFL(grid, D, dt)` which throws `std::invalid_argument` if the von Neumann stability bound is violated:
    $D \cdot dt \cdot \left(\frac{1}{h_x^2} + \frac{1}{h_y^2} + \frac{1}{h_z^2}\right) > \frac{1}{2}.$
    The message includes the actual `D`, `dt`, the per-axis `h`s, the computed CFL number, and the limit `0.5`.
  - Forwards to `explicitEuler(psi, &diff::laplacian, D * dt, n_steps)`. The `D * dt` factor folds the diffusivity into the explicit-Euler step coefficient (mathematically identical to running an unscaled Laplacian RHS at a scaled time step). `dt` and `n_steps` retain their physical meaning to the caller.
  - Doxygen block notes the CFL formula, the throw behavior, and the equivalence trick.
- No source files are added; `gridcalc` stays INTERFACE.

**Commit message:** `feat(solve): add diffuse() with explicit-diffusion CFL stability check`

## Group 3 — Tests for explicitEuler + diffuse

- Add `test/diffusion_test.cpp`:
  - `DiffusionTest.TrigEigenfunctionDecaysToAnalytic` — at `N = 32`, `D = 0.1`, `dt = 0.05`, `n_steps = 100`. `ψ_0 = sin(x) sin(y) sin(z)` on `[0, 2π]³`. Analytical `ψ(T) = exp(-3 D T) ψ_0` with `T = n_steps * dt = 5.0`. Pass when relative max-norm error vs analytical is `< 5e-2` (the discrete eigenvalue `-2(1-cos h)/h²` per axis introduces leading-order `h²/12` error per axis, integrated over 100 steps).
  - `DiffusionTest.TrigEigenfunctionSecondOrderConvergence` — sweep `N ∈ {16, 32, 64}` keeping physical `T = 5.0` fixed; at each `N`, set `dt = h² / (6 D) * 0.4` (well under the CFL limit), `n_steps = ceil(T / dt)`. Log-log slope of relative max-norm error vs `h` in `[1.8, 2.5]` (loose upper to accept the bonus from explicit-Euler's `O(dt)` time error riding on `dt ~ h²`).
  - `DiffusionTest.GaussianSanity` — ψ_0 = exp(-r²/(2σ²)) with σ=0.5 centered at (π, π, π), N=64, D=0.05, dt picked under CFL, n_steps=100. After integration, asserts:
    - All field values are finite (no NaN, no Inf).
    - Total mass (`∫ψ dV` via Phase 3 `integrate`) is conserved within `1e-10` relative.
    - The peak max-norm has decreased (Gaussian has spread out as expected).
  - `DiffusionTest.CFLViolationThrows` — call `diffuse` with `D, dt, h` chosen so that `D dt (3/h²) > 0.5`; expect `std::invalid_argument` thrown.
  - `DiffusionTest.GenericExplicitEulerOnZeroRhs` — `explicitEuler(psi, [](const Field<double>& f){ Field<double> z(f.getGrid()); return z; }, 0.1, 100)` leaves `psi` unchanged (zero RHS). Sanity check on the generic integrator.
  - `DiffusionTest.RejectsNegativeNSteps` — both `explicitEuler` and `diffuse` throw on `n_steps < 0`.
- Wire `diffusion_test.cpp` into `test/CMakeLists.txt`.

**Commit message:** `test(solve): add explicitEuler + diffuse tests (trig + Gaussian + CFL)`

## Group 4 — Phase 5 doc-notes

- `docs/user-guide/notes/phase-5-explicit-euler.md` — user-facing summary: signatures of `explicitEuler` and `diffuse`, the CFL formula and throw behavior, a worked example diffusing a trig eigenfunction and comparing to the analytical decay.
- `docs/developer-note/notes/phase-5-explicit-euler.md` — five-section structure per `specs/CLAUDE.md` step 4d:
  1. **Theory**: heat equation, eigenfunction decay on a periodic domain, von Neumann stability analysis for explicit Euler.
  2. **Math derivation**: explicit-Euler local truncation error `O(dt²)`, global error `O(dt)`; on `dt ~ h²/D` (CFL-limited) the combined space + time error is `O(h²)`; derivation of the CFL bound `D dt sum(1/h_a²) ≤ 1/2` via the Laplacian's per-axis discrete eigenvalue.
  3. **Algorithm**: per-step `tmp = rhs(psi); psi += dt * tmp`; the `D * dt` fold inside `diffuse`; throw-on-violation policy; no internal scratch persistence at Phase 5.
  4. **Design decisions**: in-place mutation rather than return-fresh; throw rather than warn; trig + Gaussian test coverage rather than periodic Gaussian image-charge sum.
  5. **References**: Strikwerda, *Finite Difference Schemes and Partial Differential Equations* (textbook with the von Neumann analysis); LeVeque (already cited from Phase 1 — additional context but not the unique citation); a permanent-URL reference for the heat-equation eigenfunction decay.

**Commit message:** `docs(notes): add Phase 5 User Guide and Developer Note notes`

## Group 5 — Bookkeeping

- Bump `CMakeLists.txt` `project(gridcalc VERSION 0.4.0)` → `0.5.0`.
- Update `test/version_test.cpp` (`VersionTest.ReturnsZeroFourZero` → `ReturnsZeroFiveZero`).
- `CHANGELOG.md`: new `## 0.5.0 — Phase 5: explicit Euler diffusion solver` entry.
- `specs/STATUS.md`: bump "Last updated"; move Phase 5 row to Done; update "Where the project stands" + "Next action" (which becomes Phase 6 — RK4 + generic time integrator interface).

**Commit message:** `chore: bump version to 0.5.0 and refresh STATUS for Phase 5`
