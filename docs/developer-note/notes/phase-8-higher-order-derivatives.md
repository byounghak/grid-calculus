# Developer Note — Phase 8: Higher-order spatial derivatives (3rd, 4th)

> Phase metadata. Spec dir:
> `specs/2026-05-02-phase-8-higher-order-derivatives/`. Adds
> `stencil::ThirdDerivative<Order>`, `stencil::FourthDerivative<Order>`
> (orders 2 and 4 each), the `diff::detail::mixedDerivative` tensor-product
> helper, 31 named partial-derivative functions across `diff/MixedPartial.hpp`
> / `diff/ThirdOrder.hpp` / `diff/FourthOrder.hpp`, and `diff::biharmonic`
> with the `diff::d4` alias. Version target: `0.8.0`.

## 1. Theory

### Multi-index partials on a Cartesian-orthogonal periodic grid

Let $\boldsymbol\alpha = (n_x, n_y, n_z) \in \mathbb{Z}_{\geq 0}^3$ be a
multi-index with total order $|\boldsymbol\alpha| = n_x + n_y + n_z$. The
partial derivative

$$
\partial^{\boldsymbol\alpha}\psi \;=\; \frac{\partial^{|\boldsymbol\alpha|}\psi}{\partial x^{n_x}\,\partial y^{n_y}\,\partial z^{n_z}}
$$

is well-defined for $\psi \in C^{|\boldsymbol\alpha|}$. On the **Cartesian-
orthogonal** grid that Phase 1 ships ($\Omega = [0,L_x) \times [0,L_y) \times
[0,L_z)$ with periodic policy and per-axis cell sizes $h_x, h_y, h_z$),
the derivative operator is **separable**:

$$
\partial^{\boldsymbol\alpha} \;=\; \partial_x^{n_x}\circ\partial_y^{n_y}\circ\partial_z^{n_z},
$$

i.e., the multi-axis derivative is the composition (in any order) of the
per-axis 1D derivatives. Discretely this means the **stencil for
$\partial^{\boldsymbol\alpha}$ on the grid factors as a tensor product**
of three 1D stencils:

$$
(D^{\boldsymbol\alpha}\psi)_{ijk} \;=\; \frac{1}{h_x^{n_x} h_y^{n_y} h_z^{n_z}}\,\sum_{a, b, c} w^{(n_x)}_a\,w^{(n_y)}_b\,w^{(n_z)}_c\,\psi_{i+a,\,j+b,\,k+c},
$$

where $\{w^{(n)}_m\}$ is the 1D central-difference weight table for
$\partial^n / \partial x^n$ at the chosen accuracy order. **The truncation
order of the multi-axis stencil equals the minimum of the per-axis
truncation orders** — proof sketch: the leading error of the per-axis 1D
stencil is $O(h^p)$ with $p$ the per-axis order; tensoring gives error
that is the sum of three such per-axis terms, each $O(h_a^{p})$, hence
the total order is $\min_a p_a$. Since this phase ships every per-axis
table at the **same** target order (`Order` parameter), the multi-axis
stencil inherits exactly that order.

### The biharmonic as a contracted scalar operator

The biharmonic is defined by **contraction** rather than by the
tensor-product form above:

$$
\nabla^4 \;=\; (\nabla^2)^2 \;=\; \biggl(\sum_a \partial_a^2\biggr)\biggl(\sum_b \partial_b^2\biggr) \;=\; \sum_a \partial_a^4 \;+\; 2\sum_{a < b} \partial_a^2\,\partial_b^2.
$$

In 3D this is **3 single-axis 4th derivatives + 3 mixed (2,2) cross-terms**
(the cross-terms picking up a factor of 2 from $\partial_a^2\partial_b^2 +
\partial_b^2\partial_a^2$ in the expansion). Phase 8's `biharmonic`
implements exactly this six-term sum directly (not via
`laplacian∘laplacian`); the "Algorithm" section explains why.

For $\psi(\mathbf{x}) = \sin(\mathbf{k}\cdot\mathbf{x})$, every
$\partial_a^2\psi = -k_a^2\psi$, so $\nabla^2\psi = -|\mathbf{k}|^2\psi$ and
$\nabla^4\psi = |\mathbf{k}|^4\psi$. This is the roadmap acceptance bar.

### Symmetry constraints on the new 1D tables

For an **odd** total derivative (ranks 1, 3, 5, …), the weight table is
**anti-symmetric**: $w_{-m} = -w_m$, so $w_0 = 0$. For an **even** total
derivative (ranks 2, 4, 6, …), the table is **symmetric**: $w_{-m} = +w_m$,
$w_0$ free. Doubling the half-width $r$ buys two more matched Taylor
orders by symmetry — odd-order error terms cancel automatically. So
central differences come in **even** truncation orders only (2, 4, 6, …)
and there is no order-3 symmetric stencil for these derivatives.

## 2. Math derivation

### Order-2 third-derivative stencil (5-point, anti-symmetric)

We seek $\partial^3 f / \partial x^3$ at $x_i$ from

$$
\partial_x^3 f(x_i) \;\approx\; \frac{1}{h^3}\sum_{m=-2}^{2} w_m\,f(x_{i+m}),
$$

with anti-symmetry $w_{-m} = -w_m$ giving $w_0 = 0$ and two free parameters
$w_1, w_2$. Taylor-expanding $f(x_i + m h)$ and grouping by derivative
order, the moment

$$
M_p \;=\; \sum_m w_m\,m^p
$$

vanishes for **even** $p$ by anti-symmetry, leaving only odd-$p$ moment
conditions. For an order-2 approximation of $\partial^3 f/\partial x^3$:

| moment | condition                  | reason                               |
|--------|----------------------------|--------------------------------------|
| $M_1 = 2 w_1 + 4 w_2$        | $= 0$  | kill $f'(x_i)$ leak                 |
| $M_3 = 2 w_1 + 16 w_2$       | $= 6 = 3!$ | match $h^3 \cdot f'''/3!$    |

Solving: $w_1 - w_2 \cdot 6 = -3$ from $M_3 - 3 M_1$, simultaneous with
$M_1 = 0$ gives $w_2 = 1/2,\ w_1 = -1$. Anti-symmetry then fills in
$w_{-1} = +1,\ w_{-2} = -1/2$. The leading uncancelled moment is $M_5$,
which sets the truncation error at $O(h^2)$.

$$
\boxed{\,\mathtt{ThirdDerivative}\langle 2\rangle\mathtt{::weights} \;=\; \{-1/2,\ +1,\ 0,\ -1,\ +1/2\}.\,}
$$

### Order-4 third-derivative stencil (7-point, anti-symmetric)

The half-width grows to $r = 3$; three free parameters $w_1, w_2, w_3$.
Order-4 accuracy demands $M_1 = 0$, $M_3 = 6$, $M_5 = 0$ (the next
uncancelled odd moment is suppressed). Solving the resulting $3 \times 3$
linear system yields

$$
\boxed{\,\mathtt{ThirdDerivative}\langle 4\rangle\mathtt{::weights} \;=\; \{1/8,\ -1,\ 13/8,\ 0,\ -13/8,\ +1,\ -1/8\}.\,}
$$

### Order-2 fourth-derivative stencil (5-point, symmetric)

For $\partial^4 f/\partial x^4$ with symmetry $w_{-m} = w_m$, the moment
conditions for an order-2 approximation are $M_0 = 0$, $M_2 = 0$,
$M_4 = 24 = 4!$. Solving gives

$$
\boxed{\,\mathtt{FourthDerivative}\langle 2\rangle\mathtt{::weights} \;=\; \{+1,\ -4,\ +6,\ -4,\ +1\}.\,}
$$

(This is the well-known 5-point fourth-difference stencil; same family
that appears as the 1D part of the 13-point biharmonic stencil and as
the second iterate of the 3-point Laplacian stencil.)

### Order-4 fourth-derivative stencil (7-point, symmetric)

Half-width $r = 3$; free $w_0, w_1, w_2, w_3$. Conditions: $M_0 = 0$,
$M_2 = 0$, $M_4 = 24$, $M_6 = 0$. Solving:

$$
\boxed{\,\mathtt{FourthDerivative}\langle 4\rangle\mathtt{::weights} \;=\; \{-1/6,\ +2,\ -13/2,\ +28/3,\ -13/2,\ +2,\ -1/6\}.\,}
$$

The truncation error is set by the next uncancelled even moment $M_8$, at
order $O(h^4)$.

### Closed-form analytical reference for the convergence sweep

For $\psi(\mathbf{x}) = \sin x \sin y \sin z$ on the periodic $[0, 2\pi)^3$
box, the partial derivative

$$
\partial^{\boldsymbol\alpha}\psi(x, y, z) \;=\; f_{n_x}(x)\,f_{n_y}(y)\,f_{n_z}(z),
$$

where $f_n$ is the $n$-th derivative of $\sin$:

$$
f_n(\xi) \;=\; \begin{cases}\sin\xi & n \equiv 0\pmod 4,\\ \cos\xi & n \equiv 1\pmod 4,\\ -\sin\xi & n \equiv 2\pmod 4,\\ -\cos\xi & n \equiv 3\pmod 4.\end{cases}
$$

This closed form drives the analytical reference in
[`test/phase_8_test.cpp`](../../../test/phase_8_test.cpp), eliminating
the need to call any of the operators under test as their own reference.
For each (operator, Order) pair, the relative max-norm error against this
reference is computed at $N \in \{16, 32, 64\}$ and the log-log slope vs
$h$ must lie in $[\text{Order} - 0.5, \text{Order} + 0.5]$.

## 3. Algorithm

### `detail::mixedDerivative<Nx, Ny, Nz, Order>` (the helper)

Defined in
[`include/gridcalc/diff/detail/MultiIndexDerivative.hpp`](../../../include/gridcalc/diff/detail/MultiIndexDerivative.hpp).
Implementation outline:

1. The internal trait `AxisStencil<N, Order>` maps the per-axis derivative
   count `N` to the corresponding 1D `gridcalc::stencil::*` table (`N = 0`
   is identity; `N ∈ {1,2,3,4}` aliases `FirstDerivative`, `Coefficients`,
   `ThirdDerivative`, `FourthDerivative` via inheritance).
2. The per-grid-point work is a triple-nested loop over the outer-product
   stencil offsets, accumulating a single double `sum`. Each loop iteration
   reads one neighbor via `psi(i + ox, j + oy, k + oz)`, which delegates
   periodic wrap to `Field::Policy`. Per-axis scaling
   $1/h_a^{N_a}$ is computed once outside the per-point loop.
3. Per-point work is $\prod_a (2 r_a + 1)$ reads, where $r_a$ is the
   radius of axis-$a$'s table (radius `0` for $N_a = 0$, radii from the
   chosen `Order`'s tables otherwise).

The named functions in `MixedPartial.hpp`, `ThirdOrder.hpp`, and
`FourthOrder.hpp` are 2-line wrappers around this helper, instantiated
with the right multi-index. There is exactly one site of correctness for
the multi-axis stencil math.

### `biharmonic<Order>(psi)` (the contracted operator)

Defined in
[`include/gridcalc/diff/Biharmonic.hpp`](../../../include/gridcalc/diff/Biharmonic.hpp).
Implementation:

1. **Single fused pass over the grid** — not `laplacian(laplacian(psi))`.
   The composition through an intermediate `Field` would inherit the
   stride-2 leading-error behavior already documented in `STATUS.md` for
   `divergence(gradient(psi))` (same convergence order but a different
   leading constant). The direct stencil sidesteps that drift and lands
   the roadmap acceptance bar at the published constant.
2. **Six-term sum** at each grid point: three single-axis 4th derivatives
   ($\partial_x^4, \partial_y^4, \partial_z^4$) using
   `stencil::FourthDerivative<Order>`, plus three mixed (2, 2) cross-terms
   ($\partial_x^2\partial_y^2, \partial_x^2\partial_z^2,
   \partial_y^2\partial_z^2$) using the outer product of two
   `stencil::Coefficients<Order>` tables along the relevant axis pair.
   Per-axis scaling: $1/h_a^4$ for the single-axis terms, $1/(h_a^2 h_b^2)$
   for the cross-terms. Per-point reads: order 2 → 27 (3 × 5 single-axis
   + 3 × 3 × 3 cross-terms with overlap); order 4 → 21 + 75 = 96 (rough).
3. **`d4<Order>(psi)`** is a 1-line alias forwarding to `biharmonic<Order>`.

### Test surface

[`test/phase_8_test.cpp`](../../../test/phase_8_test.cpp) holds:

- **Weight-table audits** — `EXPECT_DOUBLE_EQ` on every weight in
  `ThirdDerivative<2>`, `ThirdDerivative<4>`, `FourthDerivative<2>`,
  `FourthDerivative<4>` against the values derived in §2.
- **62 convergence sweeps** generated by the
  `EXPECT_PARTIAL_SLOPE_RANGE(NAME, NX, NY, NZ, ORDER)` macro — one
  per (named function, Order) pair across the 31 partials, exercising
  the user-facing API explicitly so the `<Order>` template binding is
  part of the test.
- **2 biharmonic convergence sweeps** + **4 explicit
  $\nabla^4\psi = |\mathbf{k}|^4\psi$ acceptance tests** at $N = 64$.
- **2 alias-parity tests** that verify `d4<Order>(psi)` returns the same
  data as `biharmonic<Order>(psi)` element-wise.

## 4. Design decisions

These distill the choices recorded in
`specs/2026-05-02-phase-8-higher-order-derivatives/requirements.md`.

- **Scope expansion to all 31 unique partials.** The roadmap originally
  scoped Phase 8 to "rank-2 mixed partials + biharmonic." The user
  expanded scope during the spec round to also ship single-axis 3rd & 4th
  partials and all mixed 3rd / 4th partials, on the argument that the
  multi-index machinery is the same code regardless of which subset is
  exposed and shipping the full set now avoids a follow-up phase. Roadmap
  was edited in the same branch to match.

- **`d<rank>d<axis-with-exponent>...` naming convention with lexicographic
  axis order.** The alternative — a unified `derivative<Rank, Axes...>`
  template — was rejected because at the call site the named form is
  more self-documenting and the function-list explosion is bounded
  (31 functions) and one-time. The lexicographic-ordering rule
  collapses operator-equivalent permutations (`d3dx2dy` and `d3dydx2`
  represent the same operator on smooth fields; only the lexicographic
  form is named) and keeps the public surface deterministic.

- **Direct stencils for biharmonic, not `laplacian∘laplacian`.** The
  composition would converge at the right order but with a different
  leading-error constant — the same stride-2 phenomenon already
  documented for `divergence(gradient(psi))` in Phase 2's STATUS entry.
  The roadmap acceptance bar ($\nabla^4\psi$ recovers $|\mathbf{k}|^4$
  on $\sin(\mathbf{k}\cdot\mathbf{x})$) is cleanest with the direct
  six-term stencil.

- **Sibling templates `ThirdDerivative<Order>` and `FourthDerivative<Order>`,
  not a unified `Derivative<Rank, Order>`.** Same call as Phase 2 made
  for `FirstDerivative<Order>` (sibling to Phase 1's `Coefficients`).
  Renaming the existing tables would be a public-API break; adding new
  siblings is additive. The unification meta-decision was considered for
  the second time in Phase 8 and again deferred — without a use case
  where `Rank` would be a meta-template parameter, the redesign cost
  exceeds its benefit.

- **Hard-coded weight tables, no constexpr Fornberg generator.** Phase 8
  ships eight new specializations totaling ~30 LOC of weights. Fornberg
  becomes worth its weight when a phase needs orders $\geq 6$ or many
  ranks at once; that decision moves to Phase 11 (up-to-4th-order
  gradient) at the earliest. Audit of the hard-coded values is delegated
  to `EXPECT_DOUBLE_EQ` against published Fornberg-derived references in
  the test suite.

- **Convergence-sweep grid sizes `N ∈ {16, 32, 64}`.** Same triplet as
  Phases 1, 2, 7 — keeps the slope window familiar and CI duration
  manageable. The order-4 rank-(2,2) cross stencil has radius 3 per
  axis; at $N = 16$ the halo-to-grid ratio is $3/16 = 19\%$, safe under
  `IndexPolicy::Periodic`.

## 5. References

External:

- Fornberg, B. (1988). *Generation of Finite Difference Formulas on
  Arbitrarily Spaced Grids*. **Mathematics of Computation**, 51(184),
  699–706. <https://doi.org/10.1090/S0025-5718-1988-0935077-0>. The
  generic algorithm for deriving central-difference weights of arbitrary
  derivative rank and accuracy order; used as the *theoretical* derivation
  tool for the four new tables in §2 even though the project hard-codes
  the resulting values.
- LeVeque, R. J. (2007). *Finite Difference Methods for Ordinary and
  Partial Differential Equations: Steady-State and Time-Dependent
  Problems*. SIAM. ISBN 978-0-89871-629-0.
  <https://doi.org/10.1137/1.9780898717839>. Chapter 1, Section 1.3
  ("Higher order accuracy") gives the standard 5-point fourth-difference
  stencil and the truncation analysis that yields the order-2
  fourth-derivative weights $\{1, -4, 6, -4, 1\}$.
- NIST Digital Library of Mathematical Functions, §1.10
  ("Differential Operators"). <https://dlmf.nist.gov/1.10>. Permanent
  reference for the formal definition of $\nabla^4$ as $(\nabla^2)^2$
  and the resulting six-term expansion in 3D used in
  `Biharmonic.hpp`'s implementation.

Internal cross-references (do not satisfy the "external" minimum):

- [`docs/developer-note/notes/phase-7-higher-order-stencils.md`](phase-7-higher-order-stencils.md)
  — the `Coefficients<4>` and `FirstDerivative<4>` weight tables
  derived in Phase 7 are reused unchanged by Phase 8's
  `mixedDerivative` helper for any `(N, Order)` with `N ∈ {1, 2}`.
- [`docs/developer-note/notes/phase-2-gradient-divergence.md`](phase-2-gradient-divergence.md)
  — for the documented stride-2 leading-error behavior of
  `divergence(gradient(psi))` that motivates the direct-stencil
  biharmonic implementation in §3.
