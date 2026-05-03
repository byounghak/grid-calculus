# User Guide — Phase 8: Higher-order spatial derivatives (3rd, 4th)

> Spec dir: `specs/2026-05-02-phase-8-higher-order-derivatives/`. Version: `0.8.0`.

## What this gives you

Every distinct scalar partial derivative on the periodic Cartesian grid up to
total rank 4, plus the biharmonic operator $\nabla^4$. Each operator is
templated on accuracy `Order = 2` (default) or `Order = 4`, mirroring
Phase 7's switching mechanism on `laplacian`, `gradient`, `divergence`.

```cpp
diff::dxx(psi);              // ∂²ψ/∂x² (rank 2, order 2 default)
diff::dxy<4>(psi);           // ∂²ψ/∂x∂y (rank 2, order 4)
diff::d3dx2dy(psi);          // ∂³ψ/∂x²∂y (rank 3)
diff::d3dxdydz<4>(psi);      // ∂³ψ/∂x∂y∂z (rank 3, order 4)
diff::d4dx2dy2(psi);         // ∂⁴ψ/∂x²∂y² (rank 4)
diff::d4dx2dydz<4>(psi);     // ∂⁴ψ/∂x²∂y∂z (rank 4, order 4)
diff::biharmonic(psi);       // ∇⁴ψ (rank 4, contracted)
diff::d4(psi);               // alias for biharmonic — same result
```

## Naming convention

`d<rank>d<axis-with-exponent>...`. The rules:

1. **Total rank** appears once after the leading `d` (e.g., `d3` for a
   rank-3 partial).
2. **Each distinct axis** is listed in **lexicographic order** (`x` → `y`
   → `z`), with its **exponent** (the multiplicity of differentiation along
   that axis). Exponent of `1` is omitted.
3. Concatenated as one identifier; no separators, lowercase.

| Math notation                                | Function name |
|----------------------------------------------|---------------|
| $\partial^2 \psi / \partial x^2$             | `dxx`         |
| $\partial^2 \psi / \partial y^2$             | `dyy`         |
| $\partial^2 \psi / \partial z^2$             | `dzz`         |
| $\partial^2 \psi / \partial x \partial y$    | `dxy`         |
| $\partial^2 \psi / \partial x \partial z$    | `dxz`         |
| $\partial^2 \psi / \partial y \partial z$    | `dyz`         |
| $\partial^3 \psi / \partial x^3$             | `d3dx3`       |
| $\partial^3 \psi / \partial x^2 \partial y$  | `d3dx2dy`     |
| $\partial^3 \psi / \partial x \partial y^2$  | `d3dxdy2`     |
| $\partial^3 \psi / \partial x \partial y \partial z$ | `d3dxdydz` |
| $\partial^4 \psi / \partial x^4$             | `d4dx4`       |
| $\partial^4 \psi / \partial x^3 \partial y$  | `d4dx3dy`     |
| $\partial^4 \psi / \partial x^2 \partial y^2$| `d4dx2dy2`    |
| $\partial^4 \psi / \partial x^2 \partial y \partial z$ | `d4dx2dydz` |
| $\nabla^4 \psi$ (contracted)                 | `biharmonic` *or* `d4` |

A single shorthand like `d4` (no axis suffix) denotes the **contracted
scalar operator** at that rank — $\nabla^4$ for `d4`. Phase 8 ships only the
even-rank shorthand `d4` because $\nabla^2 = $ Laplacian already has the
familiar name `laplacian` and isn't being renamed.

## API

Headers under [`include/gridcalc/diff/`](../../../include/gridcalc/diff/):

| Header                 | Operators                                                                         |
|------------------------|-----------------------------------------------------------------------------------|
| `MixedPartial.hpp`     | `dxx`, `dyy`, `dzz`, `dxy`, `dxz`, `dyz` (6 rank-2 partials)                      |
| `ThirdOrder.hpp`       | 10 rank-3 partials (`d3dx3`, …, `d3dxdydz`)                                       |
| `FourthOrder.hpp`      | 15 rank-4 partials (`d4dx4`, …, `d4dxdydz2`)                                      |
| `Biharmonic.hpp`       | `biharmonic`, `d4` alias                                                          |

Every function has the same signature shape:

```cpp
template <int Order = 2>
inline core::Field<double> NAME(const core::Field<double>& psi);
```

`Order` must be a supported even integer (`2` or `4` at Phase 8). Calling
e.g. `d3dx2dy<3>(psi)` triggers a **compile error** at the call site (the
underlying primary template is intentionally undefined), not a silent
zero-stencil bug.

### Stencil radius and per-point cost

| Operator family          | Order 2 radius / axis | Order 4 radius / axis | Per-point reads (worst case at rank 4) |
|--------------------------|-----------------------|-----------------------|----------------------------------------|
| Single-axis (`d3dx3`, `d4dx4`)     | 2 / 0 / 0       | 3 / 0 / 0             | `5` (order 2) → `7` (order 4)          |
| Two-axis (`d3dx2dy`, `d4dx2dy2`)   | 2 / 1 / 0 → 2 / 2 / 0 | 3 / 2 / 0 → 3 / 3 / 0 | up to `25` (order 2) → `49` (order 4)  |
| Three-axis (`d3dxdydz`, `d4dx2dydz`) | various       | various              | up to `27` (order 2) → `147` (order 4) |
| `biharmonic`             | 6 stencils combined   | 6 stencils combined   | order 2: ~`27`; order 4: ~`81`         |

Per-axis $1/h_a^{N_a}$ scaling means **anisotropic per-axis spacing** is
handled by the existing `Field`'s `Vec3d {hx, hy, hz}` without any
weight-table modification.

## Worked example: $\nabla^4$ on a sinusoidal eigenfunction

For $\psi(\mathbf{x}) = \sin(\mathbf{k}\cdot\mathbf{x})$ on the periodic
$[0, 2\pi)^3$ box, the biharmonic recovers $|\mathbf{k}|^4$:

$$\nabla^4 \psi = |\mathbf{k}|^4\,\psi.$$

Numerically:

```cpp
#include <cmath>
#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/diff/Biharmonic.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
using gridcalc::diff::biharmonic;

constexpr double kPi = 3.14159265358979323846;

const int N = 64;
const double L = 2.0 * kPi;
const double h = L / N;

Grid grid(N, N, N, Vec3d(h, h, h));
Field<double> psi(grid, [](double x, double y, double z) {
  return std::sin(x + 2.0 * y + z);                    // k = (1, 2, 1), |k|² = 6
});

Field<double> bh_o2 = biharmonic<2>(psi);              // ~6% relative error at N=64
Field<double> bh_o4 = biharmonic<4>(psi);              // ~3e-3 relative error at N=64

// bh ≈ |k|⁴ * psi = 36 * psi
```

For the same accuracy, order 4 lets you stay on a coarser grid and saves
substantial memory in 3D (each grid refinement is 8× the storage).

## Combining with Phase 1–7 operators

The Phase 8 operators are first-class members of `gridcalc::diff` and
compose with everything from Phases 1–7 by ordinary function composition.
Examples:

- A **biharmonic regularizer** on a free-energy functional written through
  Phase 4's `func::evaluate`: pass `[](const Field& f) { return biharmonic<4>(f); }`
  as the integrand.
- A **fourth-order PDE driver** through Phase 6's `solve::integrate`:
  build the RHS callable using `diff::biharmonic<4>` plus any other
  Phase 8 partials.

## What this does *not* do (yet)

- **No order ≥ 6 stencils.** Same call as Phase 7 — needs a Fornberg
  generator the project has not yet built.
- **No order-3 / asymmetric stencils.** Symmetric central differences have
  even truncation orders only.
- **No vector- or tensor-valued partial-derivative entry points.** The
  Phase 8 family is `Field<double> → Field<double>` only.
  Up-to-4th-order *gradient* (Phase 11) and vector / tensor fields
  (Phase 13/14) own those.
- **No spectral verification.** Phase 9 adds the FFT cross-check on every
  Phase 8 operator.
- **No non-orthogonal Bravais grids.** Phase 15 — cross-derivative terms
  appear and the per-axis-independent stencil structure breaks down.
- **No runtime-order or runtime-multi-index dispatch.** Compile-time only;
  add a thin wrapper layer downstream if config-file-driven dispatch is
  needed.

For the Taylor-matching derivation of the 3rd- and 4th-derivative weight
tables, the outer-product reasoning behind the multi-axis stencils, and
the design rationale for the d-prefix naming convention, see
[`docs/developer-note/notes/phase-8-higher-order-derivatives.md`](../../developer-note/notes/phase-8-higher-order-derivatives.md).
