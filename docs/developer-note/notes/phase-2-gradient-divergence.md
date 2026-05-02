# Developer Note — Phase 2: Gradient and divergence (scalar)

> Phase metadata. Spec dir:
> `specs/2026-05-02-phase-2-gradient-divergence/`. Implements
> `diff::gradient`, `diff::divergence`, `stencil::FirstDerivative<2>`,
> `core/EigenAliases.hpp`. Version target: `0.2.0`.
>
> *This Developer Note was written retroactively after Phase 4 (see
> specs/STATUS.md "Decisions worth knowing").*

## 1. Theory

Phase 2 adds the **first-order vector operators** to the periodic
Cartesian setting from Phase 1. For a scalar field
$\psi : \Omega \to \mathbb{R}$:

- **Gradient:** $\nabla\psi = (\partial_x \psi,\, \partial_y \psi,\, \partial_z \psi)$,
  a vector field.
- **Divergence:** $\nabla \cdot \mathbf{V} = \partial_x V_x + \partial_y V_y + \partial_z V_z$,
  a scalar field.

Both decompose into independent 1D first-derivative operators along each
Cartesian axis, mirroring the Cartesian decomposition Phase 1 used for
the Laplacian. The key new ingredient is a **first-derivative
central-difference stencil**, sibling to the second-derivative table
that already existed.

A small but consequential observation: composing first-derivatives gives
an effective **stride-2 second-derivative stencil**, *not* the 3-point
stencil that `laplacian` uses. Concretely,

$$
(D_1 \circ D_1)f(x_i) \;=\; \frac{D_1 f(x_{i+1}) - D_1 f(x_{i-1})}{2h} \;=\; \frac{f(x_{i+2}) - 2 f(x_i) + f(x_{i-2})}{4 h^2},
$$

while

$$
D_2 f(x_i) \;=\; \frac{f(x_{i-1}) - 2 f(x_i) + f(x_{i+1})}{h^2}.
$$

Both converge to $\partial_x^2 f$ as $h \to 0$, but at finite $h$ they
have **different leading-order constants**. The Phase 2 round-trip test
`DivergenceTest.RoundTripWithGradientMatchesLaplacian` therefore
compares `divergence(gradient(ψ))` to the *analytical* Laplacian, not
to `laplacian(ψ)`.

## 2. Math derivation

### 1D first-derivative central-difference stencil

Subtracting Taylor expansions:

$$
f(x + h) - f(x - h) \;=\; 2 h f'(x) + \tfrac{h^3}{3} f'''(x) + O(h^5).
$$

Dividing by $2h$:

$$
\partial_x f(x) \;\approx\; \frac{f(x + h) - f(x - h)}{2 h}
\;=\; f'(x) + \tfrac{h^2}{6} f'''(x) + O(h^4),
$$

with truncation error $\tfrac{h^2}{6} f'''(x)$. The 3-element table
`stencil::FirstDerivative<2>::{offsets, weights}` encodes weights
$\{-\tfrac{1}{2},\, 0,\, +\tfrac{1}{2}\}$ unscaled; the consumer
multiplies by $1/h$. Phase 7 adds the `Order = 4` specialization
alongside `Coefficients<4>`.

### Discrete eigenvalue on a single trig

For $\psi(x) = \sin(k x)$ sampled on the periodic grid:

$$
\frac{\sin(k(x + h)) - \sin(k(x - h))}{2 h} \;=\; \frac{2 \cos(k x)\,\sin(k h)}{2 h} \;=\; \frac{\sin(k h)}{h} \cos(k x).
$$

So the discrete first-derivative operator has eigenvalue
$\tilde{k} = \sin(k h)/h$ on the trig basis function, vs the continuous
$k$. Taylor: $\sin(k h)/h = k - k^3 h^2 / 6 + O(h^4)$, so the **relative
error** is $-k^2 h^2 / 6 + O(h^4)$. At $N = 32$, $h = 2\pi/32 \approx 0.196$,
$k = 1$: relative error $\approx -0.64\%$. This is the bound the
`GradientTest.RecoversTrigGradient` test relaxes to `1.5e-2` for
headroom.

### 3D gradient and divergence on Cartesian coordinates

Apply `FirstDerivative<2>` along each axis independently:

$$
(\nabla_h \psi)_{ijk} \;=\; \begin{pmatrix}
\frac{\psi_{i+1, j, k} - \psi_{i-1, j, k}}{2 h_x} \\[4pt]
\frac{\psi_{i, j+1, k} - \psi_{i, j-1, k}}{2 h_y} \\[4pt]
\frac{\psi_{i, j, k+1} - \psi_{i, j, k-1}}{2 h_z}
\end{pmatrix},
\qquad
(\nabla_h \cdot \mathbf{V})_{ijk} \;=\; \sum_{a \in \{x,y,z\}} \frac{V_{a,\,\text{shift}(+\hat{a})} - V_{a,\,\text{shift}(-\hat{a})}}{2 h_a}.
$$

Each per-axis term is independent; each inherits the 1D truncation error
$\tfrac{h_a^2}{6} \partial_a^3 \psi$ derived above.

### Composition: $\nabla_h \cdot \nabla_h \psi$ vs $\nabla^2_h \psi$

Composing the gradient and divergence stencils yields a **stride-2
3-point stencil** along each axis:

$$
\bigl[(\nabla_h \cdot \nabla_h)\psi\bigr]_{ijk} \;=\; \sum_a \frac{\psi_{ijk + 2 \hat{a}} - 2\psi_{ijk} + \psi_{ijk - 2 \hat{a}}}{(2 h_a)^2}.
$$

Per axis, the discrete eigenvalue on $\sin(k_a x_a)$ is

$$
-\frac{2(1 - \cos(2 k_a h_a))}{(2 h_a)^2} \;=\; -\frac{\sin^2(k_a h_a)}{h_a^2} \;=\; -k_a^2 + \frac{k_a^4 h_a^2}{3} + O(h_a^4),
$$

while the direct 3-point Laplacian of Phase 1 has eigenvalue (per axis)

$$
-\frac{2(1 - \cos(k_a h_a))}{h_a^2} \;=\; -k_a^2 + \frac{k_a^4 h_a^2}{12} + O(h_a^4).
$$

Both are $O(h^2)$ and converge to the analytical eigenvalue $-k_a^2$,
but the **leading constants differ** — the composed operator's error is
4× larger than the direct Laplacian's. Hence the round-trip test
compares to the analytical eigenvalue: `divergence(gradient(ψ))`
converges to $\nabla^2 \psi$, but at finite $h$ it does *not* equal
`laplacian(ψ)`.

## 3. Algorithm

### `FirstDerivative<2>` table

[`include/gridcalc/stencil/FirstDerivative.hpp`](../../../include/gridcalc/stencil/FirstDerivative.hpp)
follows the same primary-undefined + `Order = 2` specialization pattern
as Phase 1's `Coefficients<Order>`. Weights are `{-0.5, 0, 0.5}` and the
`weights[1] = 0` is included so the loop body can be the same generic
"sum weights × shifted samples" structure as the Laplacian's.

### `gradient`

[`include/gridcalc/diff/Gradient.hpp`](../../../include/gridcalc/diff/Gradient.hpp)
allocates a fresh `Field<Vec3d>` of the same shape and loops
i-fastest. At each point, it accumulates three independent
`FirstDerivative<2>`-weighted sums (one per axis) into three doubles,
then packs them as `Vec3d(dx*inv_hx, dy*inv_hy, dz*inv_hz)`. Periodic
wrap is delegated to the input `Field<double>`'s `IndexPolicy::Periodic`
via `operator()`.

### `divergence`

[`include/gridcalc/diff/Divergence.hpp`](../../../include/gridcalc/diff/Divergence.hpp)
allocates a fresh `Field<double>` and at each point reads only the
matching component from the shifted `Field<Vec3d>`: `V(i+s, j, k)(0)`
for the x-derivative, `V(i, j+s, k)(1)` for the y-derivative,
`V(i, j, k+s)(2)` for the z-derivative. The cross-component reads
(e.g., `V(i+s, j, k)(1)`) that a generic Vec3d-stencil would touch
are intentionally absent — divergence is component-along-axis only.

Algorithmic complexity for both is $O(N_x N_y N_z)$. Allocation cost:
gradient creates a `Field<Vec3d>` of $24 N_x N_y N_z$ bytes (3 doubles
per point) plus a fresh `Field<double>` for the result of any
follow-on divergence; divergence creates a `Field<double>` for the
output ($8 N_x N_y N_z$ bytes).

## 4. Design decisions

- **Gradient returns `Field<Vec3d>`, not three `Field<double>`s.**
  Single field-of-vectors matches the roadmap. `Field<T>` is templated
  from Phase 1, so this is a one-line return — no wrapper struct, no
  helper. Per-component decomposition can be added later as an additive
  helper without changing the canonical API. See
  [Phase 2 requirements.md](../../../specs/2026-05-02-phase-2-gradient-divergence/requirements.md).
- **`Vec3d` lives in `core/EigenAliases.hpp`.** Phase 2 lifted the alias
  out of `core/Grid.hpp` (where Phase 1 declared it) once a second
  operator needed it. The fully-qualified name `gridcalc::core::Vec3d`
  is unchanged. Future tensor / matrix typedefs (`Mat3d`, `Tensor3` in
  Phase 13/14) will land in this same header.
- **Divergence accepts a single `Field<Vec3d>`.** Symmetric with the
  gradient return; no per-component overload at Phase 2. A 3-component
  overload would just duplicate the implementation; if a downstream
  caller needs it, it can be added then.
- **`stencil::FirstDerivative<Order>` is a sibling template, not part
  of `Coefficients<Order>`.** Phase 1's `Coefficients<Order>` is
  documented as the **second-derivative** table; reusing the name for
  first-derivative coefficients would silently change semantics for
  Phase 1 callers and break the public contract. A unifying redesign
  (e.g., `Coefficients<Derivative, Order>`) is deferred to Phase 7
  when both tables gain `Order = 4` specializations.
- **Round-trip identity test compares to the analytical Laplacian, not
  the discrete one.** As §2 derives, $D_1 \circ D_1 \ne D_2$ at finite
  $h$. This is a feature of the Phase 2 design that is easy to misread
  as a bug; the test name `RoundTripWithGradientMatchesLaplacian` is
  about the continuum limit, not the discrete one.

## 5. References

[1] LeVeque, R. J. (2007). *Finite Difference Methods for Ordinary and
Partial Differential Equations*. SIAM. ISBN 978-0-89871-783-9.
<https://doi.org/10.1137/1.9780898717839>. §1.2 covers first-derivative
central differences and the truncation analysis used in §2; §1.5
discusses the multi-stencil compositions that motivate the
$D_1 \circ D_1 \ne D_2$ analysis.

[2] Fornberg, B. (1988). "Generation of Finite Difference Formulas on
Arbitrarily Spaced Grids." *Mathematics of Computation*, 51(184),
699–706. <https://doi.org/10.1090/S0025-5718-1988-0935077-0>. The
canonical reference for deriving central-difference weights at
arbitrary order; the `Order = 4` extension in Phase 7 will use this
algorithm directly.

[3] Trefethen, L. N. (2000). *Spectral Methods in MATLAB*. SIAM.
<https://doi.org/10.1137/1.9780898719598>. §3 on the relationship
between discrete first-derivatives and the trigonometric eigenvalue
$\sin(k h)/h$ used in §2.

[4] Eigen 3 documentation: `Eigen::Vector3d` and uninitialized
default-construction.
<https://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html>
(permanent URL). Documents that `Eigen::Matrix` default-constructed
fixed-size objects leave their coefficients uninitialized in release
builds; this is why Phase 2's tests broadcast-init `Field<Vec3d>` with
`Vec3d::Zero()` rather than relying on the `Field(grid)` constructor's
documented "zero-initialized" promise (which holds for arithmetic `T`
but not for Eigen types).
