# Developer Note — Phase 1: Periodic 3D scalar grid + Laplacian

> Phase metadata. Spec dir: `specs/2026-05-01-phase-1-laplacian/`.
> Implements `core::Grid`, `core::Field<T, Policy = Periodic>`,
> `core::IndexPolicy::Periodic`, `stencil::Coefficients<2>`, and
> `diff::laplacian`. Version target: `0.1.0`.
>
> *This Developer Note was written retroactively after Phase 4, when the
> Phase-3-onward doc-notes rule was extended to cover Phases 1 and 2 as
> well (see specs/STATUS.md "Decisions worth knowing").*

## 1. Theory

Phase 1 stands up the **discrete approximation of $\nabla^2$** on a
periodic Cartesian domain $\Omega = [0, L_x) \times [0, L_y) \times [0, L_z)$
with $L_a = N_a h_a$ for $a \in \{x, y, z\}$. A scalar field $\psi : \Omega \to \mathbb{R}$ is sampled at the lattice points
$\mathbf{x}_{ijk} = (i h_x,\, j h_y,\, k h_z)$ for $0 \le i < N_x$,
$0 \le j < N_y$, $0 \le k < N_z$.

The continuous Laplacian on this domain has tensor-product trig
eigenfunctions: for $\psi_{\mathbf{k}}(\mathbf{x}) = \prod_a \sin(k_a x_a)$
with integer wavenumbers $k_a \in \mathbb{Z}_{>0}$,

$$
\nabla^2 \psi_{\mathbf{k}} \;=\; -(k_x^2 + k_y^2 + k_z^2)\,\psi_{\mathbf{k}}.
$$

Eigenfunctions of the discrete Laplacian are the **same trig products**
sampled at the lattice; the discrete eigenvalues differ from the
continuous ones by a quantifiable $O(h^2)$ correction (derived in §2).
The 3D Laplacian decomposes as a sum of 1D second derivatives along each
axis (Cartesian-orthogonal ⇒ no mixed partials), so the implementation
applies the 1D 2nd-order central-difference stencil along each axis
independently and sums.

## 2. Math derivation

### 1D second-derivative central-difference stencil

Taylor-expanding $f(x \pm h)$:

$$
f(x + h) = f(x) + h f'(x) + \tfrac{h^2}{2} f''(x) + \tfrac{h^3}{6} f'''(x) + \tfrac{h^4}{24} f''''(x) + O(h^5),
$$

$$
f(x - h) = f(x) - h f'(x) + \tfrac{h^2}{2} f''(x) - \tfrac{h^3}{6} f'''(x) + \tfrac{h^4}{24} f''''(x) + O(h^5).
$$

Adding cancels the odd derivatives:

$$
f(x + h) + f(x - h) - 2 f(x) \;=\; h^2 f''(x) + \tfrac{h^4}{12} f''''(x) + O(h^6).
$$

Dividing by $h^2$ gives the **2nd-order accurate stencil** that
`stencil::Coefficients<2>` encodes:

$$
\partial^2_x f(x) \;\approx\; \frac{f(x - h) - 2 f(x) + f(x + h)}{h^2}
\;=\; f''(x) + \tfrac{h^2}{12} f''''(x) + O(h^4),
$$

with **truncation error $-\tfrac{h^2}{12} f''''(x)$** (using the
discrete-minus-continuous sign convention). Weights `{1, -2, 1}` divided
by $h^2$, encoded as the 3-element table
`Coefficients<2>::{offsets, weights}` per
[`include/gridcalc/stencil/CentralDifference.hpp`](../../../include/gridcalc/stencil/CentralDifference.hpp).
The `Order = 4` specialization is deferred to Phase 7.

### 3D Laplacian via Cartesian decomposition

The 3D operator is a sum of three 1D operators, one per axis, each
scaled by its own $1/h_a^2$:

$$
\nabla^2_h \psi_{ijk} \;=\; \sum_{a \in \{x, y, z\}} \frac{1}{h_a^2} \sum_{m=-1}^{+1} w_m\, \psi(\mathbf{x}_{ijk} + m h_a \hat{\mathbf{e}}_a),
$$

with $w = (1, -2, 1)$. This is implemented as three nested loops in
[`include/gridcalc/diff/Laplacian.hpp`](../../../include/gridcalc/diff/Laplacian.hpp)
that each accumulate into a per-point scalar before writing to the
output `Field`.

### Discrete eigenvalue on a trig product

For $\psi_{ijk} = \sin(k_x x_i) \sin(k_y y_j) \sin(k_z z_k)$ on the
periodic grid (with each $k_a$ a positive integer so that the function
is grid-periodic):

$$
\frac{\sin(k_a (x_a - h_a)) - 2\sin(k_a x_a) + \sin(k_a (x_a + h_a))}{h_a^2} \;=\; -\frac{2(1 - \cos(k_a h_a))}{h_a^2}\,\sin(k_a x_a),
$$

so $\psi_{\mathbf{k}}$ is an exact eigenfunction of $\nabla^2_h$ with
eigenvalue

$$
\lambda_h(\mathbf{k}) \;=\; -\sum_a \frac{2(1 - \cos(k_a h_a))}{h_a^2}.
$$

Taylor-expanding $1 - \cos(k_a h_a) = (k_a h_a)^2/2 - (k_a h_a)^4/24 + O(h^6)$ gives

$$
\lambda_h(\mathbf{k}) \;=\; -\sum_a k_a^2 \;+\; \frac{1}{12} \sum_a k_a^4 h_a^2 \;+\; O(h^4),
$$

i.e. the **continuous eigenvalue $-\sum_a k_a^2$ plus an $O(h^2)$
correction** with leading coefficient $\tfrac{1}{12} \sum_a k_a^4$. For
the test `LaplacianTest.RecoversTrigEigenvalue` at $N = 32$, $h = 2\pi/32$,
$\mathbf{k} = (1, 2, 1)$:

$$
\frac{|\lambda_h - \lambda|}{|\lambda|} \;\approx\; \frac{(1^4 + 2^4 + 1^4)\,h^2}{12 \cdot (1^2 + 2^2 + 1^2)} \;=\; \frac{18 h^2}{72} \;=\; \frac{h^2}{4} \;\approx\; 9.6 \times 10^{-3}.
$$

This matches the test's tolerance of `1.5e-2` with healthy headroom and
is the per-axis-independent error model the convergence-sweep test
exploits at $N \in \{16, 32, 64\}$.

## 3. Algorithm

### Grid

[`core::Grid`](../../../include/gridcalc/core/Grid.hpp) is a small POD
holding $(N_x, N_y, N_z)$ and `Vec3d cell_size` $(h_x, h_y, h_z)$.
`getLinearIndex(i, j, k)` returns the i-fastest offset `i + Nx*(j + Ny*k)`.
Construction validates that every extent and every cell-size component
is strictly positive.

### Field<T, Policy>

[`core::Field<T, Policy>`](../../../include/gridcalc/core/Field.hpp)
holds a `Grid` by value and a `std::vector<T>` of length
$N_x N_y N_z$ in i-fastest layout. Three constructors:

1. `Field(grid)` — zero-initialized via `std::vector<T>(N)` value-init.
2. `Field(grid, value)` — broadcast.
3. `Field(grid, f)` — SFINAE-enabled when `f` is invocable as
   `f(double, double, double) -> T`; samples `f` at the Cartesian
   coordinates `(i*hx, j*hy, k*hz)`.

`operator()(i, j, k)` always passes its arguments through
`Policy::wrap` before calling `Grid::getLinearIndex`. For the default
`IndexPolicy::Periodic`, this is the modular wrap
`((i % N) + N) % N`. The wrap is always applied, never special-cased
on whether `i` is in range — branchless modular arithmetic in the hot
path is the right tradeoff because it eliminates the branch predictor
miss that a range check would incur on every stencil access.

### Laplacian operator

[`diff::laplacian`](../../../include/gridcalc/diff/Laplacian.hpp) loops
over all grid points (z-slowest, x-fastest, matching the storage
layout). At each point, it accumulates three independent
`Coefficients<2>`-weighted sums — one per axis — and returns
`dxx*inv_hx2 + dyy*inv_hy2 + dzz*inv_hz2`. Because `operator()` already
wraps via the periodic policy, no separate boundary handling code path
is needed; the same loop body works for boundary and interior points.

Algorithmic complexity is $O(N_x N_y N_z)$. The output is a fresh
`Field<double>`, allocated once at the start.

## 4. Design decisions

- **`Field<T, Policy>` is templated on the policy from Phase 1.** Only
  `Periodic` exists at Phase 1; `Dirichlet` and `Neumann` arrive in
  Phase 18. The policy template parameter was the *only* piece of
  forward-looking design admitted at Phase 1 — specifically so Phase 18
  is an additive change, not a public-API break. See
  [Phase 1 requirements.md](../../../specs/2026-05-01-phase-1-laplacian/requirements.md).
- **i-fastest storage layout.** Linear offset
  `i + Nx*(j + Ny*k)`. Matches PocketFFT / FFTW (Phase 9) and Eigen's
  column-major default. Public-contract status from Phase 1 forward;
  Phase 20 benchmarks would be invalidated by a layout flip later.
- **`stencil::Coefficients<Order>` is a template; primary undefined.**
  Only `Order = 2` is specialized at Phase 1; instantiating an
  unsupported `Order` trips a compile error at the call site rather
  than silently producing zero output. Phase 7 adds the `Order = 4`
  specialization with no call-site refactor.
- **Eigen propagated as a `SYSTEM` include.** Required so
  `-Wconversion` (and the rest of the strict warning set) stays active
  on the project's own code without firing on Eigen's NEON typecasting
  headers. `WARNINGS_AS_ERRORS=ON` in CI breaks without this. See
  the `tech-stack.md` Eigen section.
- **`Vec3d` originally lived in `core/Grid.hpp`.** Phase 2 relocated
  the alias to `core/EigenAliases.hpp` once `gradient` and `divergence`
  also needed it; the fully-qualified name `gridcalc::core::Vec3d` is
  unchanged.

## 5. References

[1] LeVeque, R. J. (2007). *Finite Difference Methods for Ordinary and
Partial Differential Equations*. SIAM. ISBN 978-0-89871-783-9.
<https://doi.org/10.1137/1.9780898717839>. §1.2 on second-derivative
central differences and the truncation analysis used in §2 of this note;
§2 on Cartesian-product extension of 1D stencils to higher dimensions.

[2] Trefethen, L. N. (2000). *Spectral Methods in MATLAB*. SIAM. ISBN
978-0-89871-465-4. <https://doi.org/10.1137/1.9780898719598>. Chapter 1
covers the spectral relationship between trig basis functions and finite
differences on periodic domains, the framing used in §1 of this note.

[3] cppreference: `std::vector` value-initialization semantics.
<https://en.cppreference.com/w/cpp/container/vector/vector> (permanent
URL). Confirms the zero-initialization guarantee for `Field(grid)` when
`T` is `double`. Note: this guarantee does *not* hold for Eigen
fixed-size types (the default constructor leaves coefficients
uninitialized) — Phase 2's tests work around this by explicitly
broadcast-initializing `Field<Vec3d>` with `Vec3d::Zero()`.
