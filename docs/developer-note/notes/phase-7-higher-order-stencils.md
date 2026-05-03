# Developer Note — Phase 7: Higher-order accuracy stencils (4th-order)

> Phase metadata. Spec dir:
> `specs/2026-05-03-phase-7-higher-order-stencils/`. Adds
> `stencil::Coefficients<4>`, `stencil::FirstDerivative<4>`, and
> templates `diff::laplacian`, `diff::gradient`, `diff::divergence`
> on an `Order` parameter (default `2`). Version target: `0.7.0`.

## 1. Theory

Central-difference stencils approximate a derivative at $x_i$ as a
weighted sum of samples at $x_{i-r}, x_{i-r+1}, \ldots, x_{i+r}$. The
weights $w_m$ are chosen so the weighted sum matches the target
derivative through some order in the Taylor expansion of $f$ around
$x_i$.

Two **symmetry constraints** organize the family:

- For an **odd** derivative ($\partial^k$ with $k$ odd), the stencil is
  **anti-symmetric**: $w_{-m} = -w_m$, so the center weight $w_0$ is
  exactly zero. This is why `FirstDerivative<2>::weights[1] = 0` and
  `FirstDerivative<4>::weights[2] = 0`.
- For an **even** derivative, the stencil is **symmetric**:
  $w_{-m} = +w_m$, with no constraint on $w_0$.

Doubling the half-width $r$ buys two more matched Taylor orders by
symmetry — odd-order error terms cancel automatically. So central
differences come in **even** truncation orders only: 2, 4, 6, ...
There is no order-3 symmetric central-difference stencil for these
derivatives; an order-3 method would have to be asymmetric (upwind),
which is a different design space.

For the Phase 1 `Grid` (Cartesian-orthogonal with per-axis cell sizes
$h_x, h_y, h_z$), the 3D Laplacian decomposes into independent
per-axis 1D second-derivative stencils:

$$
\nabla^2_h \psi_{ijk} \;=\; \sum_{a \in \{x,y,z\}} \frac{1}{h_a^2}\,\sum_{m=-r}^{+r} w_m\,\psi(\mathbf{x}_{ijk} + m\,h_a\,\hat{\mathbf{e}}_a),
$$

with the same weight table $\{w_m\}$ on every axis (because the weights
are unitless — they only depend on the integer offsets $\{-r, \ldots, +r\}$,
not on $h_a$). This is what lets one `Coefficients<Order>` table drive
the operator on **anisotropic Cartesian-orthogonal** grids without any
per-axis customization.

## 2. Math derivation

### Order-4 second-derivative stencil

We want $\partial^2 f / \partial x^2$ at $x_i$ from the symmetric stencil

$$
\partial^2_x f(x_i) \;\approx\; \frac{1}{h^2}\,\sum_{m=-2}^{2} w_m\,f(x_{i+m})
\;=\; \frac{1}{h^2}\bigl[w_{-2} f_{i-2} + w_{-1} f_{i-1} + w_0 f_i + w_1 f_{i+1} + w_2 f_{i+2}\bigr].
$$

Symmetry: $w_{-m} = w_m$. Three free parameters: $w_0$, $w_1 = w_{-1}$,
$w_2 = w_{-2}$.

Taylor-expand around $x_i$:

$$
f(x_i + m\,h) \;=\; \sum_{p \ge 0} \frac{(m\,h)^p}{p!}\,f^{(p)}(x_i).
$$

The weighted sum is

$$
\sum_m w_m\,f(x_i + m\,h) \;=\; \sum_{p \ge 0} \frac{h^p}{p!}\,f^{(p)}(x_i)\,\underbrace{\sum_m w_m\,m^p}_{\,M_p\,}.
$$

By symmetry $M_p = 0$ for **odd** $p$ — the odd Taylor terms cancel
automatically. The remaining moment conditions for an order-4
approximation of $\partial^2 f / \partial x^2$ are:

| moment | required value | reason                                         |
|--------|----------------|------------------------------------------------|
| $M_0 = w_0 + 2 w_1 + 2 w_2$  | `0`  | kill the constant (no $f$ leak)               |
| $M_2 = 2 w_1 + 8 w_2$       | `2`  | match $h^2 \cdot f''(x_i) \cdot (1/2)$ → coefficient `2` after multiplying by `2!/h²` |
| $M_4 = 2 w_1 + 32 w_2$      | `0`  | kill the $h^4$ term to lift accuracy from $O(h^2)$ to $O(h^4)$ |

Solving:

- From $M_2$ and $M_4$: subtract → $24 w_2 = -2$, so $w_2 = -1/12$.
- Back into $M_2$: $2 w_1 + 8 \cdot (-1/12) = 2$ → $w_1 = (2 + 2/3)/2 = 4/3$.
- From $M_0$: $w_0 = -2(w_1 + w_2) = -2(4/3 - 1/12) = -2 \cdot 5/4 \cdot ... = -5/2$. (Direct: $w_0 = -2 \cdot 4/3 - 2 \cdot (-1/12) = -8/3 + 1/6 = -16/6 + 1/6 = -15/6 = -5/2$.)

Hence

$$
\boxed{\;w = \bigl(-\tfrac{1}{12},\ \tfrac{4}{3},\ -\tfrac{5}{2},\ \tfrac{4}{3},\ -\tfrac{1}{12}\bigr)\;}
$$

with truncation error from the next non-zero even moment
$M_6 = 2 w_1 + 128 w_2 = 8/3 - 32/3 = -8$, giving leading error
$-\tfrac{h^4}{6!}\cdot M_6 \cdot f^{(6)}(x_i) = +\tfrac{8 h^4}{720}\,f^{(6)}(x_i) = \tfrac{h^4}{90}\,f^{(6)}(x_i)$.

So the local truncation error is $\tfrac{h^4}{90}\,f^{(6)}(x_i) = O(h^4)$,
which is the **order-4 accuracy** advertised. Compare with the order-2
stencil's $\tfrac{h^2}{12}\,f^{(4)}$.

### Order-4 first-derivative stencil

Same procedure for $\partial f / \partial x$ at $x_i$:

$$
\partial_x f(x_i) \;\approx\; \frac{1}{h}\,\sum_{m=-2}^{2} w_m\,f(x_{i+m}),
$$

with **anti-symmetry** $w_{-m} = -w_m$ (so $w_0 = 0$, $w_{-1} = -w_1$,
$w_{-2} = -w_2$). The non-trivial moment conditions are now odd $p$:

| moment | required value | reason                                       |
|--------|----------------|----------------------------------------------|
| $M_1 = 2 w_1 + 4 w_2$       | `1`  | match $h \cdot f'(x_i)$ → after $1/h$ scale |
| $M_3 = 2 w_1 + 16 w_2$      | `0`  | kill the $h^3$ term to lift accuracy        |

Solving:

- Subtracting $M_3 - M_1$: $12 w_2 = -1$ → $w_2 = -1/12$.
- Back into $M_1$: $2 w_1 + 4 \cdot (-1/12) = 1$ → $w_1 = (1 + 1/3)/2 = 2/3$.

Hence

$$
\boxed{\;w = \bigl(\tfrac{1}{12},\ -\tfrac{2}{3},\ 0,\ \tfrac{2}{3},\ -\tfrac{1}{12}\bigr)\;}
$$

Truncation error from $M_5$:
$M_5 = 2 w_1 + 64 w_2 = 4/3 - 16/3 = -4$, giving leading error
$-\tfrac{h^4}{5!}\cdot M_5 \cdot f^{(5)}(x_i) = \tfrac{4 h^4}{120}\,f^{(5)}(x_i) = \tfrac{h^4}{30}\,f^{(5)}(x_i) = O(h^4)$.

Compare with the order-2 first-derivative stencil's $\tfrac{h^2}{6}\,f^{(3)}$.

### Why the per-axis-independent structure works on anisotropic grids

The two derivations above produce **unitless weight tables** —
$\{-1/12, 4/3, -5/2, 4/3, -1/12\}$ for $\partial^2/\partial x^2$ and
$\{1/12, -2/3, 0, 2/3, -1/12\}$ for $\partial/\partial x$. Each table
encodes the linear combination of *index-shifted samples* that
approximates the derivative; the spacing $h$ enters only as a scalar
postmultiplication ($1/h^2$ or $1/h$).

For 3D on a Cartesian-orthogonal grid, the Laplacian decomposes axis-by-axis
because the operator $\nabla^2 = \partial_x^2 + \partial_y^2 + \partial_z^2$
has no cross-derivative terms — the basis vectors are mutually
perpendicular. So we apply the same 1D weight table along each axis
independently, scaled by that axis's own cell size. The tables don't
need to be regenerated per axis.

The order-2 Phase 1 Laplacian already used this structure; Phase 7
just substitutes a longer table.

## 3. Algorithm

The implementation is the smallest change consistent with the design:

### Stencil tables

[`stencil/CentralDifference.hpp`](../../../include/gridcalc/stencil/CentralDifference.hpp)
adds the explicit `Coefficients<4>` specialization with `radius = 2`,
`offsets = {-2, -1, 0, 1, 2}`, and the weights derived above.
[`stencil/FirstDerivative.hpp`](../../../include/gridcalc/stencil/FirstDerivative.hpp)
adds the analogous `FirstDerivative<4>` specialization. Both follow
the existing primary-undefined pattern: `Coefficients<3>` /
`FirstDerivative<3>` etc. are still undefined and trip a compile error
at the call site.

### Operator templating

The three Phase 1/2 operator headers gain a `template <int Order = 2>`
parameter on the function. The body changes by exactly one line — the
hardcoded stencil type becomes the parameterized one:

```diff
- using Coeffs = stencil::Coefficients<2>;
+ using Coeffs = stencil::Coefficients<Order>;
```

(and analogously for `FirstDerivative<Order>` in `Gradient.hpp` and
`Divergence.hpp`). The per-axis loop already iterates `m` over
`Coeffs::offsets.size()`, so a 5-element table just runs the loop one
more iteration per axis without any structural change.

### Backward compatibility

C++17 resolves a function-template call with all template parameters
defaulted using the no-template-arg syntax: `laplacian(field)` calls
`laplacian<2>(field)` because `Order` defaults to `2`. So every
Phase 1/2 caller (and every Phase 4–6 transitive caller, like
`func::evaluate(... f(_, _, lap_psi))` and `solve::diffuse`) continues
to compile and produce bit-identical output. The Phase 1/2 test files
are untouched; the Phase 7 test file (`test/higher_order_test.cpp`)
exercises the new `<4>` overload.

### Per-point cost

Order 4 reads `2 * radius + 1 = 5` samples per axis instead of `3`.
Per-grid-point cost is therefore `~5/3 ≈ 1.67×` the order-2 cost —
in exchange for the truncation error dropping from $O(h^2)$ to $O(h^4)$,
which at the typical test resolution `N = 32` corresponds to
**10–100× absolute-accuracy improvement** (verified by
`HigherOrderTest.*MoreAccurateThanOrder2`). For most applications this
is a win; the Phase 20 performance pass will quantify exact slowdown
on the benchmark suite.

## 4. Design decisions

- **Compile-time `<Order>` over runtime parameter.** Matches the
  project's existing pattern (`stencil::Coefficients<Order>` is already
  templated this way; Phase 3's `Pairwise`/`Kahan`, Phase 4's arity
  dispatch, Phase 6's `ExplicitEuler`/`RK4`). No runtime branch in the
  inner loop. Compile-time error on unsupported orders.
- **Hard-coded weights, derivation in the Developer Note.** Two
  specializations, ~10 LOC each. Deferred constexpr Fornberg until
  Phase 8/11 actually require many more orders. See `requirements.md`
  for the full rationale.
- **Default `Order = 2` rather than removing the default.** Backward
  compatibility — Phase 1/2/4/5/6 callers don't need to change. C++17
  function-template default-argument semantics make this clean.
- **Apply to all three operators in one phase.** Matches the
  roadmap deliverable list. Divergence is templated for symmetry even
  though the roadmap text only names Laplacian and gradient — the
  operators are tightly coupled (divergence reads gradient's
  `FirstDerivative` table) and skipping divergence would force a
  Phase 7.5 follow-up.
- **No `Order` parameter on `solve::diffuse`.** `solve::diffuse` is the
  thin convenience driver for $D \nabla^2 \psi$; it currently calls
  `diff::laplacian(...)` without a template arg, which now resolves to
  order 2. A future phase can promote the diffusion driver to take an
  `Order` parameter if the use case justifies it; today the workaround
  (call `solve::integrate` directly with a custom RHS that uses
  `laplacian<4>`) is documented in the User Guide.

## 5. References

[1] Fornberg, B. (1988). "Generation of Finite Difference Formulas on
Arbitrarily Spaced Grids." *Mathematics of Computation*, 51(184),
699–706. <https://doi.org/10.1090/S0025-5718-1988-0935077-0>. The
canonical reference for deriving central-difference weights at
arbitrary order on arbitrary stencils. Used here as the *theoretical*
derivation framework even though the two specific specializations are
hard-coded; the derivation in §2 above is a hand-evaluation of
Fornberg's algorithm at order 4 with stencil offsets `{-2, -1, 0, 1, 2}`.

[2] LeVeque, R. J. (2007). *Finite Difference Methods for Ordinary and
Partial Differential Equations*. SIAM. ISBN 978-0-89871-783-9.
<https://doi.org/10.1137/1.9780898717839>. §1.5 covers higher-order
finite differences and lists the same `{-1/12, 4/3, -5/2, 4/3, -1/12}`
stencil for $\partial^2/\partial x^2$ in Table 1.1.

[3] DLMF (Digital Library of Mathematical Functions), §3.4
"Differentiation."
<https://dlmf.nist.gov/3.4> (permanent URL). NIST's online reference
includes the central-difference weights for both first and second
derivatives at orders 2, 4, 6, 8 in tabular form; matches the
derivation here.

[4] Trefethen, L. N. (2000). *Spectral Methods in MATLAB*. SIAM.
<https://doi.org/10.1137/1.9780898719598>. Chapter 1 motivates the
choice between increasing stencil order vs refining the grid (the
classic FD-vs-spectral tradeoff); the Phase 7 User Guide note's "When
to pick order 4" guidance is in the spirit of this discussion.
