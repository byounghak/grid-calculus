# Developer Note — Phase 6: RK4 + generic time integrator interface

> Phase metadata. Spec dir:
> `specs/2026-05-02-phase-6-rk4-integrator/`. Adds
> `solve::RK4` integrator, `solve::Integrator.hpp` aggregator,
> tag-dispatched `solve::integrate(...)`. Refactors `solve::diffuse`
> to template on the integrator tag and `solve::explicitEuler` to be
> a thin wrapper. Version target: `0.6.0`.

## 1. Theory

The Runge–Kutta family generalizes forward Euler by evaluating the RHS
at carefully chosen intermediate points within each step and combining
the results to match more terms of the Taylor expansion of the exact
flow. Higher-order RK schemes have:

- **Higher accuracy** per step (local truncation error `O(dt^{p+1})`,
  global error `O(dt^p)`, where `p` is the method order).
- **Larger stability regions** in the complex plane, so the method
  tolerates larger `dt` for stiff-but-not-too-stiff problems.
- **Multiple RHS evaluations** per step (one per stage), making them
  more expensive per step.

The classic 4-stage Runge–Kutta method (RK4) has order `p = 4` and uses
4 RHS evaluations per step. Its Butcher tableau is

$$
\begin{array}{c|cccc}
0   &       &       &       &       \\
1/2 & 1/2   &       &       &       \\
1/2 & 0     & 1/2   &       &       \\
1   & 0     & 0     & 1     &       \\
\hline
    & 1/6   & 1/3   & 1/3   & 1/6
\end{array}
$$

Per-step update for $\dot\psi = \mathrm{rhs}(\psi)$:

$$
\begin{aligned}
k_1 &= \mathrm{rhs}(\psi^n) \\
k_2 &= \mathrm{rhs}(\psi^n + (dt/2)\,k_1) \\
k_3 &= \mathrm{rhs}(\psi^n + (dt/2)\,k_2) \\
k_4 &= \mathrm{rhs}(\psi^n + dt\,k_3) \\
\psi^{n+1} &= \psi^n + (dt/6)(k_1 + 2 k_2 + 2 k_3 + k_4).
\end{aligned}
$$

For the heat equation $\partial_t\psi = D\nabla^2\psi$, RK4 advances each
trig eigenfunction by the truncated Taylor series of
$e^{D\,\lambda_h\,dt}$ — exact through fourth order in the per-step
exponent.

## 2. Math derivation

### RK4's order-4 accuracy

Apply RK4 to a linear scalar test problem $\dot y = z\,y$ with
$z = D\,\lambda_h$ (the discrete eigenvalue along one trig mode).
Stage values evaluate to:

$$
\begin{aligned}
k_1 &= z\,y, \\
k_2 &= z\,(y + (dt/2)\,k_1) = z\,y\,(1 + z\,dt/2), \\
k_3 &= z\,(y + (dt/2)\,k_2) = z\,y\,(1 + z\,dt/2 + (z\,dt)^2/4), \\
k_4 &= z\,(y + dt\,k_3) = z\,y\,(1 + z\,dt + (z\,dt)^2/2 + (z\,dt)^3/4).
\end{aligned}
$$

Substituting into the update:

$$
y^{n+1} \;=\; y^n\,\Bigl[1 + z\,dt + \tfrac{(z\,dt)^2}{2} + \tfrac{(z\,dt)^3}{6} + \tfrac{(z\,dt)^4}{24}\Bigr].
$$

This is exactly the **fourth-order Taylor truncation** of
$e^{z\,dt}$. The next term, $(z\,dt)^5/120$, is the local truncation
error — `O(dt^5)` per step, hence `O(dt^4)` global error after `T/dt`
steps.

### Stability region intersection with the negative real axis

The amplification factor $g(z\,dt)$ above must satisfy $|g| \le 1$ for
stability. For real negative $z\,dt = -x$ ($x > 0$),

$$
g(-x) \;=\; 1 - x + \tfrac{x^2}{2} - \tfrac{x^3}{6} + \tfrac{x^4}{24}.
$$

Setting $g(-x) = -1$ and solving numerically yields $x \approx 2.7853$.
For $x \in [0,\,2.7853]$, $g(-x) \in [-1,\,1]$ and the method is
stable; beyond that, $|g| > 1$ and the discrete solution amplifies.

For the heat equation with $z = D\,\lambda_h$ and the worst-case
discrete eigenvalue $\lambda_h^{\max} = -\sum_a 4/h_a^2$, the
RK4 stability bound is

$$
4\,D\,dt\,\sum_a \frac{1}{h_a^2} \;\le\; 2.7853, \qquad
\boxed{\,D\,dt\,\sum_a \frac{1}{h_a^2} \;\le\; 0.6963.\,}
$$

This is the value encoded as `RK4::diffusionCFLLimit`. Compared with
forward Euler's `0.5`, RK4 is `~1.39×` looser — modest but useful when
the spatial grid is fine and Euler is at the edge of its CFL.

### Why combined refinement at CFL-limited dt is spatial-dominated

When `dt = c · h^2 / D` (CFL-limited), the time-error term scales
differently for the two integrators:

- **Euler** time error per step: $\sim |z|^2\,dt^2 / 2 = O(h^4)$ per
  step; over $T/dt = O(h^{-2})$ steps, total time error $O(h^2)$.
  Combined with spatial error $O(h^2)$, total error is **$O(h^2)$**.
- **RK4** time error per step: $|z|^5\,dt^5 / 120 = O(h^{10})$ per
  step; over $O(h^{-2})$ steps, total time error $O(h^8)$. Spatial
  error remains $O(h^2)$. Combined error is **$O(h^2)$** — also
  spatial-dominated, but the constant is smaller.

So a combined-refinement convergence sweep (the user-chosen test
shape) shows slope `~2` for **both** integrators; the RK4 vs Euler
benefit appears in the **constant**, not the **slope**. The
`IntegratorTest.RK4LessErrorThanEulerAtSameDt` test demonstrates that
constant differs by at least 10× when the time-error contribution is
fully isolated from the spatial-error contribution (using a pure-ODE
manufactured solution with no Laplacian).

The roadmap-mandated **dt^4 slope** is verified by
`IntegratorTest.RK4TimeAccuracyOrder4` — also a pure-ODE test, since
the slope is hidden behind the spatial floor for any PDE-based test
under CFL-limited stepping.

## 3. Algorithm

### Tag types

[`include/gridcalc/solve/Integrator.hpp`](../../../include/gridcalc/solve/Integrator.hpp)
is an aggregator that pulls in
[`ExplicitEuler.hpp`](../../../include/gridcalc/solve/ExplicitEuler.hpp)
and [`RK4.hpp`](../../../include/gridcalc/solve/RK4.hpp). Each
per-integrator header defines its own tag struct:

```cpp
struct ExplicitEuler { static constexpr double diffusionCFLLimit = 0.5;    };
struct RK4           { static constexpr double diffusionCFLLimit = 0.6963; };
```

This co-location avoids a circular include — if the tags lived in
`Integrator.hpp` and `Integrator.hpp` included the per-integrator
headers, then those per-integrator headers (which need the tag types)
couldn't include `Integrator.hpp` back without a cycle.

### Generic `integrate` overloads

Each per-integrator header defines its `integrate(...)` overload tagged
on its own tag struct. So the dispatch is plain function-overload
resolution at the call site; no `if constexpr` ladder, no virtual
dispatch.

`integrate(... RK4)` (in `RK4.hpp`) implements the four-stage update
using a small `detail::axpyFresh(a, scale, b)` helper that returns
`a + scale * b` element-wise as a fresh `Field`. Per step, this
allocates four stage outputs (`k1, k2, k3, k4`) plus three intermediate
states (`psi + (dt/2) k1`, `psi + (dt/2) k2`, `psi + dt k3`) — about
six `Field<double>` allocations per step at $N_x N_y N_z$ doubles each.
Phase 20 will swap this for a persistent scratch pool with a
non-allocating RHS callback shape.

The final accumulation is fused into one element-wise loop over
`(k1 + k4)/6 + (k2 + k3)/3`, reducing reads from the four `k_i` arrays
without an intermediate sum field.

### `explicitEuler` thin wrapper

[`ExplicitEuler.hpp`](../../../include/gridcalc/solve/ExplicitEuler.hpp)
keeps Phase 5's free function name:

```cpp
template <class Rhs>
inline void explicitEuler(core::Field<double>& psi, Rhs&& rhs,
                          double dt, int n_steps) {
  integrate(psi, std::forward<Rhs>(rhs), dt, n_steps, ExplicitEuler{});
}
```

Phase 5 callers compile unchanged; the
`DiffusionTest.GenericExplicitEulerOnZeroRhs` test continues to exercise
the wrapper directly and `IntegratorTest.GenericIntegrateDispatchEulerEqualsExplicitEuler`
verifies the wrapper produces bit-identical output to the underlying
`integrate(..., ExplicitEuler{})` call.

### `diffuse` becomes a template

[`Diffusion.hpp`](../../../include/gridcalc/solve/Diffusion.hpp):

```cpp
template <class Integrator = ExplicitEuler>
inline void diffuse(core::Field<double>& psi, double D, double dt,
                    int n_steps, Integrator tag = {}) {
  /* validate args */
  detail::checkDiffusionCFL<Integrator>(psi.getGrid(), D, dt);
  integrate(
      psi,
      [](const core::Field<double>& f) { return diff::laplacian(f); },
      D * dt,
      n_steps,
      tag);
}
```

The CFL check reads `Integrator::diffusionCFLLimit` so each integrator
gets its own bound. `typeid(Tag).name()` is included in the error
message for diagnostics.

## 4. Design decisions

- **Tag dispatch over class hierarchy.** Compile-time dispatch matches
  the rest of the project (Phase 3's `Pairwise`/`Kahan`, Phase 4's
  `evaluate` arity). No virtual-call overhead in the per-step inner
  loop. The tag's `static constexpr` properties (CFL limit) are
  available to higher-level drivers without any indirection.
- **Tag types co-located with their integrators.** Avoids the circular
  include described in §3 above. The downside — a caller writing
  `solve::ExplicitEuler{}` needs to have included
  `solve/ExplicitEuler.hpp` (or `solve/Integrator.hpp` which transitively
  pulls it) — is acceptable. The `diffuse` driver pulls
  `Integrator.hpp` so any tag works through that entry point.
- **`explicitEuler` kept as wrapper rather than removed.** Phase 5
  callers and tests continue working unchanged. The two equivalent
  forms are documented as a stylistic choice.
- **`diffuse` templated on integrator rather than two separate drivers.**
  One driver, one CFL-check implementation, parameterized on the tag.
  Adding a third integrator later (e.g., Phase 19's implicit BDF) needs
  only a new tag struct and a new `integrate(... BDF)` overload — no
  change to `diffuse`.
- **Pure-ODE convergence tests for `dt^4` slope.** PDE-based tests at
  CFL-limited `dt` are spatial-dominated for both integrators (§2);
  the user's chosen combined-refinement test shows slope `~2` (the
  expected combined order). To honor the roadmap's literal `Δt^4`
  acceptance, an additional test uses a pure linear-decay ODE
  `dpsi/dt = -α ψ` with no spatial discretization, where the time-error
  slope is unambiguous.
- **No persistent scratch pool yet.** Per-step allocation cost (~6
  fields for RK4) is documented as a known cost. Phase 20 (performance)
  is the right place to introduce a non-allocating RHS shape and a
  reusable scratch.

## 5. References

[1] Hairer, E., Nørsett, S. P., & Wanner, G. (2008). *Solving Ordinary
Differential Equations I: Nonstiff Problems* (2nd revised ed.).
Springer Series in Computational Mathematics, vol. 8. Springer-Verlag.
ISBN 978-3-540-56670-0. <https://doi.org/10.1007/978-3-540-78862-1>.
Chapter II.1 covers the classical Runge–Kutta methods including the
specific 4-stage scheme used here, with the Butcher-tableau derivation
in §II.2 and the order-condition analysis in §II.3.

[2] Butcher, J. C. (2016). *Numerical Methods for Ordinary Differential
Equations* (3rd ed.). Wiley. ISBN 978-1-119-12150-3.
<https://doi.org/10.1002/9781119121534>. The canonical modern reference
for Runge–Kutta theory; §3 gives the order-condition framework that
proves RK4's local truncation error is `O(dt^5)`.

[3] Wikipedia: *Runge–Kutta methods* — *The Runge–Kutta method* (RK4
section) and *Stability* subsection.
<https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods> (permanent
URL). Includes the explicit form of the RK4 update used in §3 and the
stability-region plot showing the `~−2.7853` real-axis intersection.

[4] LeVeque, R. J. (2007). *Finite Difference Methods for Ordinary and
Partial Differential Equations*. SIAM. ISBN 978-0-89871-783-9.
<https://doi.org/10.1137/1.9780898717839>. §7.4 covers RK methods
applied to method-of-lines discretizations of parabolic PDEs, the
framework used in `solve::diffuse`. Cited from earlier phases as well;
this chapter is new context for Phase 6.
