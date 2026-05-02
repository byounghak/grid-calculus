# Requirements: Phase 5 — Explicit Euler diffusion solver

## Roadmap reference

> ## Phase 5 — Explicit Euler diffusion solver
>
> - **Goal.** Time integration enters the picture: $\partial_t\psi = D\nabla^2\psi$.
> - **Deliverables.**
>   - `solve/ExplicitEuler.hpp` — fixed time-step integrator parameterized by an RHS callable.
>   - `solve/Diffusion.hpp` — convenience driver `diffuse(psi, D, dt, n_steps)`.
>   - Stability check: warns/errors when $D\,dt/h^2 > 0.5/d$ (CFL for $d$-D explicit diffusion).
> - **Acceptance.** Numerical solution to a Gaussian initial condition matches the analytic solution within finite-difference error after 100 steps.

## Decisions made for this feature

- **Two-layer API: `solve::explicitEuler(psi, rhs, dt, n_steps)` + `solve::diffuse(psi, D, dt, n_steps)`.** Matches the roadmap's two-deliverable layout. The integrator is generic over an RHS callable `Field<double> rhs(const Field<double>&)`; the diffusion driver supplies `diff::laplacian` as that RHS. The two are independent public APIs in the `gridcalc::solve` namespace.
- **In-place mutation: `void diffuse(Field<double>& psi, …)` and `void explicitEuler(Field<double>& psi, …)`.** Standard for iterative time integration. The caller copies the IC first if they want to keep the unchanged starting state. This diverges from the Phase 1–4 operator style (which returns fresh fields), but is the right pattern for a multi-step integrator where allocating one output field per step would inflate memory churn.
- **CFL stability check throws `std::invalid_argument`.** `solve::diffuse` validates the von Neumann bound $D \cdot dt \cdot (1/h_x^2 + 1/h_y^2 + 1/h_z^2) \le 1/2$ before integration begins and throws if violated. The error message includes the offending values and the limit. Caller catches and adjusts `dt` if they want a custom retry. Plain `solve::explicitEuler` does **not** check anything stability-related — it is generic over the RHS, so it cannot know the relevant stability bound.
- **Per-axis CFL form rather than the roadmap's `D dt / h² > 0.5/d` shorthand.** `Grid` supports anisotropic per-axis cell sizes (Phase 1). The form $D \cdot dt \cdot \sum_a (1/h_a^2) \le 1/2$ is the von Neumann derivation's actual statement and reduces to the roadmap's $D dt / h^2 \le 1/(2d)$ when $h_a = h$ for all axes. Same content, anisotropy-aware.
- **`diffuse` folds `D` into the explicit-Euler step coefficient.** Internally, `diffuse(psi, D, dt, n_steps)` calls `explicitEuler(psi, &diff::laplacian, D * dt, n_steps)`. Mathematically: $\psi^{n+1} = \psi^n + dt \cdot (D \nabla^2 \psi^n) = \psi^n + (D \cdot dt) \cdot \nabla^2 \psi^n$. This avoids allocating an intermediate scaled Field per step inside the RHS lambda. The user-visible `dt` and `n_steps` retain their physical meaning; the trick is implementation-private.
- **Test set.** Five tests:
  1. Trig eigenfunction matching analytical decay at `N = 32, T = 5.0` (relative error `< 5e-2`).
  2. Convergence sweep on the trig eigenfunction over `N ∈ {16, 32, 64}` with `dt ∝ h²` (CFL-limited); log-log slope `~2`.
  3. Gaussian sanity (no NaN, mass conserved, peak decreased) — the roadmap's literal Gaussian-IC requirement, but framed as a sanity test rather than a convergence comparison (the periodic-image-charge analytical reference is too heavy a lift for a Phase 5 deliverable).
  4. CFL-violation test: `diffuse` throws on out-of-range parameters.
  5. Generic-integrator sanity: `explicitEuler(psi, zero-RHS, dt, 100) leaves psi unchanged`.
- **Per-step allocation accepted.** Each `explicitEuler` step allocates one fresh `Field<double>` (the RHS-callable's return). For a 100-step run on a `64³` field, this is 100 alloc/free pairs of ~2 MB each — significant but not catastrophic, and Phase 20 (performance) is the right place to optimize via persistent scratch buffers + a non-allocating RHS callback shape (e.g., `void rhs(const Field&, Field& out)`). Anticipating that change at Phase 5 would expand API surface speculatively.
- **`int n_steps` rather than `std::size_t`.** Matches Phase 1's convention for grid extents (`int Nx`); makes "0 steps" representable for testing without unsigned-wrap surprises. Negative values are rejected at runtime.
- **Version bump: `0.5.0`** on merge. New public symbols (`solve::explicitEuler`, `solve::diffuse`) and new public namespace `gridcalc::solve`; no existing symbol semantics change.

## Non-goals

- **No implicit / RK4 / generic time-integrator interface.** Phase 6's deliverable. A class-based or virtual-interface integrator at Phase 5 would either be discarded or refactored when Phase 6 lands.
- **No adaptive time stepping.** Out of scope for the roadmap; would require local-error estimation and is most natural after RK4.
- **No periodic-Gaussian image-charge analytical reference.** The 3D Jacobi-theta image sum is heavy infrastructure for one test; replaced by the trig-eigenfunction analytical comparison plus a Gaussian sanity test that doesn't need an analytical answer.
- **No vector or tensor diffusion.** Phase 11+ vector-field operators arrive then.
- **No automatic CFL adjustment.** `diffuse` throws; the caller adjusts and retries. No auto-shrinking of `dt`.
- **No threading.** The accumulation loop is serial. Phase 20 parallelizes (and inherits the Phase 3 determinism contract for the `func::integrate` call inside `DiffusionTest.GaussianSanity`).

## Deferred questions

- **Per-step allocation removal via out-parameter RHS callable.** Deferred to Phase 20 (performance pass). Will require either a second callable shape `void rhs(const Field&, Field&)` or an Eigen-style expression-template fusion (Phase 11+).
- **Generic integrator interface (`solve::Integrator`).** Phase 6 designs and ships this; explicit-Euler will be refactored to consume it then.
- **Implicit-diffusion variant.** Phase 19's deliverable, after Phase 18 introduces non-periodic boundary policies.
- **Adaptive time stepping.** Not on the current roadmap; would naturally land alongside RK4 in Phase 6 if pursued.
