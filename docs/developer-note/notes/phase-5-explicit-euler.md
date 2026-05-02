# Developer Note — Phase 5: Explicit Euler diffusion solver

> Phase metadata. Spec dir:
> `specs/2026-05-02-phase-5-explicit-euler/`. Implements
> `solve::explicitEuler` (generic) and `solve::diffuse` (heat-equation
> driver with CFL check). Version target: `0.5.0`.

## 1. Theory

Phase 5 introduces **time integration** to the library. The continuous
problem is the heat equation on a periodic Cartesian domain:

$$
\partial_t \psi(\mathbf{x}, t) \;=\; D\, \nabla^2 \psi(\mathbf{x}, t),
\qquad \psi(\mathbf{x}, 0) = \psi_0(\mathbf{x}),
$$

with $\psi : \Omega \times [0, T] \to \mathbb{R}$, $D > 0$, and periodic
boundary conditions on $\Omega = [0, L_x) \times [0, L_y) \times [0, L_z)$.

**Forward (explicit) Euler** discretizes time as
$\psi^{n+1} = \psi^n + dt \cdot \mathrm{rhs}(\psi^n)$. For the heat
equation, $\mathrm{rhs}(\psi) = D \nabla^2 \psi$.

The eigenfunctions of the continuous Laplacian on this periodic domain
are tensor-product trig modes (Phase 1 §1). They are also eigenfunctions
of the **discrete** Laplacian $\nabla^2_h$, which is convenient: the
continuous solution and the semi-discrete (continuous-time) solution and
the fully-discrete (forward-Euler) solution can all be written
analytically as products of factors of the form $e^{-\lambda t}$ or
$(1 - dt\,\lambda)^n$. This makes the trig-IC test in §2 a direct
algebraic comparison.

**Stability** is the new ingredient. The discrete Laplacian is bounded
below (its smallest eigenvalue is the most negative); explicit Euler is
only conditionally stable, and the bound on $dt$ comes from a von
Neumann analysis of the discrete amplification factor.

## 2. Math derivation

### Discrete eigenvalues of the Laplacian (recap from Phase 1)

For $\psi_{\mathbf{k}} = \prod_a \sin(k_a x_a)$ with integer
wavenumbers $k_a$, the Phase-1 stencil gives (per axis):

$$
\frac{\psi_{ijk + \hat{a}} - 2 \psi_{ijk} + \psi_{ijk - \hat{a}}}{h_a^2}\Big|_{\psi_{\mathbf{k}}} \;=\; -\frac{2(1 - \cos(k_a h_a))}{h_a^2}\,\psi_{\mathbf{k}}.
$$

Summing over axes,

$$
\nabla^2_h \psi_{\mathbf{k}} \;=\; \lambda_h(\mathbf{k})\,\psi_{\mathbf{k}}, \qquad \lambda_h(\mathbf{k}) \;=\; -\sum_a \frac{2(1 - \cos(k_a h_a))}{h_a^2}.
$$

### Forward-Euler amplification factor

Because $\psi_{\mathbf{k}}$ is an eigenfunction of $\nabla^2_h$, the
discrete update on this mode is a scalar:

$$
\psi^{n+1}_{\mathbf{k}} \;=\; \psi^n_{\mathbf{k}} + dt \cdot D \lambda_h(\mathbf{k})\,\psi^n_{\mathbf{k}} \;=\; \bigl(1 + D\,dt\,\lambda_h(\mathbf{k})\bigr)\,\psi^n_{\mathbf{k}}.
$$

The **amplification factor** is $g(\mathbf{k}) = 1 + D\,dt\,\lambda_h(\mathbf{k})$.
Since $\lambda_h(\mathbf{k}) \le 0$ for all $\mathbf{k}$ (the Laplacian is
negative semi-definite), $g \le 1$ automatically. The lower bound
$g \ge -1$ is the binding stability requirement:

$$
1 + D\,dt\,\lambda_h(\mathbf{k}) \;\ge\; -1
\quad\Longleftrightarrow\quad
D\,dt\,|\lambda_h(\mathbf{k})| \;\le\; 2.
$$

Maximizing over $\mathbf{k}$ — each $1 - \cos(k_a h_a)$ is maximized at
$2$ (when $k_a h_a = \pi$) — gives the **CFL bound**:

$$
\max_{\mathbf{k}} |\lambda_h(\mathbf{k})| \;=\; \sum_a \frac{4}{h_a^2}, \qquad \boxed{\,D\,dt\,\sum_a \frac{1}{h_a^2} \;\le\; \frac{1}{2}.\,}
$$

For isotropic $h$, this reduces to $D\,dt/h^2 \le 1/(2d)$ — the form
quoted in the roadmap. The implementation uses the per-axis sum because
`Grid` allows anisotropic spacing.

### Convergence on the trig IC

For $\psi_0 = \sin x \sin y \sin z$ on $[0, 2\pi]^3$ (so $\mathbf{k} =
(1, 1, 1)$), the analytical solution is
$\psi(t) = e^{-3 D t}\,\psi_0$. After $n$ explicit-Euler steps the
discrete solution is

$$
\psi^n \;=\; (1 + D\,dt\,\lambda_h)^n\,\psi_0, \qquad \lambda_h \;=\; -\sum_a \frac{2(1 - \cos h_a)}{h_a^2}.
$$

For $h_a = h$ small, $\lambda_h = -3 + h^2/4 + O(h^4)$ (using
$1 - \cos h = h^2/2 - h^4/24 + \ldots$ per axis times 3 axes; verifying:
$3 \cdot (h^2/2 - h^4/24) \cdot 2/h^2 = 3 - h^2/4 + O(h^4)$, so
$\lambda_h = -3 + h^2/4 + O(h^4)$).
Hence $1 + D\,dt\,\lambda_h = 1 - 3 D\,dt + (D\,dt)\,h^2/4 + O(dt h^4)$.

Comparing $(1 + D\,dt\,\lambda_h)^n$ to $e^{-3 D t}$ at $t = n\,dt$:

- The **time-discretization** error from $(1 + x)^n$ vs $e^{n x}$ is
  $O(dt)$ globally for explicit Euler.
- The **spatial** error from $\lambda_h - (-3) = h^2/4 + O(h^4)$
  contributes an extra $O(h^2)$.
- When the CFL is at its limit, $dt \sim h^2/D$, so the time error
  inherits the same $O(h^2)$ scaling as the spatial error. The combined
  error is **$O(h^2)$** under fixed-$T$ refinement.

This is the basis for `DiffusionTest.TrigEigenfunctionSecondOrderConvergence`.
The log-log slope test allows up to $2.5$ rather than $2.2$ because the
$O(dt)$ time-error term contributes a small bonus when $dt$ is below the
CFL limit (we use $dt = 0.4 \cdot h^2/(6D)$, well under the bound).

## 3. Algorithm

### Generic integrator

[`include/gridcalc/solve/ExplicitEuler.hpp`](../../../include/gridcalc/solve/ExplicitEuler.hpp).

```cpp
template <class Rhs>
inline void explicitEuler(core::Field<double>& psi, Rhs&& rhs, double dt, int n_steps);
```

A SFINAE-style `static_assert` validates that `Rhs` is callable as
`Field<double>(const Field<double>&)`. Each step:

1. `Field<double> tmp = rhs(psi)` — one fresh allocation.
2. `psi[n] += dt * tmp[n]` for `n` in `[0, N)` — raw-pointer loop.

The integrator runs entirely on `psi.data()` and `tmp.data()`; it
deliberately does **not** use `Field::operator()` here because the
periodic-wrap policy contributes nothing inside the in-place
accumulation (we touch every element exactly once, in storage order).

### Diffusion driver

[`include/gridcalc/solve/Diffusion.hpp`](../../../include/gridcalc/solve/Diffusion.hpp).

```cpp
inline void diffuse(core::Field<double>& psi, double D, double dt, int n_steps);
```

Validates `D, dt, n_steps >= 0`, then calls `checkExplicitDiffusionCFL`
(which throws on violation), then forwards to
`explicitEuler(psi, &diff::laplacian, D * dt, n_steps)`. The `D * dt`
fold is the algebraic identity
$\psi^{n+1} = \psi^n + dt(D\nabla^2 \psi^n) = \psi^n + (D\,dt)\nabla^2 \psi^n$;
it lets the explicit-Euler step accept an unscaled Laplacian RHS without
a per-step element-wise scaling loop, at the cost of `dt` losing its
physical-time meaning **inside `explicitEuler`'s frame**. The user-visible
`dt` in the `diffuse` API is unaffected.

### Per-step allocation

Each `explicitEuler` step allocates one fresh `Field<double>` (the RHS
return value) and frees it at the end of the step. For a 100-step run on
a $64^3$ field, this is 100 alloc/free pairs of ~2 MB each. The
allocator's freelist makes this much cheaper than 100 fresh `mmap`s in
practice (most modern `malloc` implementations recycle the buffer).
Phase 20 will switch to a persistent scratch via a non-allocating RHS
shape (e.g., `void rhs(const Field&, Field& out)`), which is the
standard approach for tight time-stepping loops.

### CFL check

`checkExplicitDiffusionCFL(grid, D, dt)` computes
$\mathrm{cfl} = D\,dt \sum_a (1/h_a^2)$ and throws
`std::invalid_argument` if $\mathrm{cfl} > 0.5$. The message includes
the actual values for diagnostics — easy to re-compute the right `dt`
from the printed line.

## 4. Design decisions

- **In-place mutation** rather than return-fresh. Time integration
  doesn't naturally produce a single output; the standard pattern is to
  evolve a state vector. Returning a fresh `Field<double>` per call (let
  alone per step) would double memory churn for no gain. Phase 1–4
  operators return fresh because they're snapshot-style queries; Phase 5
  is the first time the iterator pattern dominates.
- **Throw on CFL violation** rather than warn. Running a 100-step explicit
  Euler past its stability bound produces NaN / Inf within a few steps;
  silent overflow is a debugging trap. A throw with a descriptive
  message lets the caller catch, halve `dt`, and retry, which is the
  expected workflow.
- **No generic integrator interface yet.** Phase 6 introduces
  `solve::Integrator` and `RK4`. Anticipating that interface here would
  ship a class hierarchy that gets refactored a phase later. Keep
  `explicitEuler` as a free function; the Phase 6 refactor will adapt
  it to consume the new interface.
- **Trig-IC convergence test rather than periodic-Gaussian.** The
  roadmap's literal "Gaussian initial condition" requires summing
  Gaussian image charges over all periodic shifts (a 3D Jacobi theta
  function) for an analytical reference; that's out of proportion to
  the rest of Phase 5's surface area. The trig eigenfunction has a
  closed-form analytical solution that's just $e^{-3 D t}$ times the IC,
  which makes the comparison trivial to read in the test. The Gaussian
  IC still appears in `GaussianSanity` as a no-NaN / mass-conserved /
  peak-decreased smoke test.
- **Per-axis CFL form**. Anisotropy-aware. See `requirements.md`.

## 5. References

[1] Strikwerda, J. C. (2004). *Finite Difference Schemes and Partial
Differential Equations* (2nd ed.). SIAM. ISBN 978-0-89871-639-9.
<https://doi.org/10.1137/1.9780898717938>. §6.1 covers the von Neumann
analysis of forward-Euler discretization of the heat equation, the
derivation reproduced in §2 of this note. §10 gives the multi-dimensional
extension and the per-axis sum form of the CFL bound.

[2] LeVeque, R. J. (2007). *Finite Difference Methods for Ordinary and
Partial Differential Equations*. SIAM. ISBN 978-0-89871-783-9.
<https://doi.org/10.1137/1.9780898717839>. §9.2 (forward-Euler heat
equation) and §10.2 (von Neumann analysis); cited from Phase 1 as well
but the heat-equation chapter is new context here.

[3] Wikipedia: *Heat equation* — *Periodic boundary conditions and
Fourier series solution*.
<https://en.wikipedia.org/wiki/Heat_equation#The_heat_equation_on_the_circle>
(permanent URL). Documents the Fourier-series solution
$\psi(t) = \sum_{\mathbf{k}} \hat\psi_{\mathbf{k}}(0)\,e^{-D|\mathbf{k}|^2 t}\,e^{i\mathbf{k}\cdot\mathbf{x}}$
and the eigenfunction decay used in the trig-IC test.

[4] Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P.
(2007). *Numerical Recipes: The Art of Scientific Computing* (3rd ed.).
Cambridge University Press. ISBN 978-0-521-88068-8. §17.2 covers the
"forward-Euler is conditionally stable" framing and the practical
recommendation to use higher-order time integrators (RK4 in Phase 6) or
implicit methods (Phase 19) for non-trivial diffusion problems.
