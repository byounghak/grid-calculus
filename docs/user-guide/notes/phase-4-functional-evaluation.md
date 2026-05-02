# User Guide — Phase 4: Functional evaluation

> Spec dir: `specs/2026-05-02-phase-4-functional-evaluation/`. Version: `0.4.0`.

## What this gives you

A primitive that evaluates a **functional** of a scalar field — an integral
of an arbitrary user-supplied integrand $f(\psi, \nabla\psi, \nabla^2\psi)$
over the periodic Cartesian domain:

$$
F[\psi] \;=\; \int_{\Omega} f\bigl(\psi(\mathbf{x}),\, \nabla\psi(\mathbf{x}),\, \nabla^2\psi(\mathbf{x})\bigr)\, dV.
$$

`evaluate` builds on top of Phase 1's `laplacian`, Phase 2's `gradient`, and
Phase 3's `integrate`: it materializes whichever derivatives your `f`
actually needs, evaluates `f` pointwise into a temporary `Field<double>`,
and reduces the result with the same pairwise / Kahan tag dispatch from
Phase 3.

## API

Header:
[`include/gridcalc/func/Functional.hpp`](../../../include/gridcalc/func/Functional.hpp).

```cpp
namespace gridcalc::func {

template <class F, class Tag = Pairwise>
double evaluate(const core::Field<double>& psi, F&& f, Tag tag = {});

}  // namespace gridcalc::func
```

Your callable `f` may take **any one of three signatures**, picked at
compile time by `std::is_invocable_r_v`:

| Signature                                         | Computes |
|---------------------------------------------------|----------|
| `f(double psi)`                                   | only `ψ` — gradient and Laplacian are **skipped** |
| `f(double psi, const Vec3d& grad)`                | `ψ` and `∇ψ` — Laplacian is **skipped** |
| `f(double psi, const Vec3d& grad, double lap)`    | full triplet `ψ`, `∇ψ`, `∇²ψ` |

Return type must be convertible to `double`. The reducer tag (`Pairwise{}`
or `Kahan{}`) is forwarded to `func::integrate`.

## Worked example: Ginzburg–Landau free energy

The Ginzburg–Landau free-energy density is

$$
f_{\mathrm{GL}}(\psi, \nabla\psi) \;=\; \tfrac{1}{2} |\nabla\psi|^2 \;+\; W (\psi^2 - 1)^2,
$$

with $W$ a positive constant. Note that $f_{\mathrm{GL}}$ does **not**
depend on $\nabla^2\psi$, so writing it with the 2-arg signature lets
`evaluate` skip the Laplacian computation.

### Setup

```cpp
#include <cmath>
#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/func/Functional.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
using gridcalc::func::evaluate;

constexpr double kPi = 3.14159265358979323846;

const int N = 64;
const double L = 2.0 * kPi;
const double h = L / N;
const double W = 1.0;

Grid grid(N, N, N, Vec3d(h, h, h));
Field<double> psi(grid, [](double x, double y, double z) {
  return std::sin(x) * std::sin(y) * std::sin(z);
});
```

### The integrand

```cpp
auto gl_density = [W](double p, const Vec3d& g) {
  const double dpsi2 = p * p - 1.0;
  return 0.5 * g.squaredNorm() + W * dpsi2 * dpsi2;
};
```

### Evaluate

```cpp
const double F = evaluate(psi, gl_density);
```

For this `ψ` with `W = 1`, the analytical reference is
$F = (3/2)\pi^3 + (411/64)\pi^3 = (507/64)\pi^3 \approx 245.62$.
At `N = 64`, `evaluate` matches this to better than 0.5% relative
error; the discrete error is dominated by the squared-gradient term and
converges at order 2 in $h$ (verified by
[`test/functional_test.cpp`](../../../test/functional_test.cpp)'s
`GinzburgLandauSecondOrderConvergence` test).

For the why-this-formula-converges-at-order-2 background — derivation,
error analysis, references — see
[`docs/developer-note/notes/phase-4-functional-evaluation.md`](../../developer-note/notes/phase-4-functional-evaluation.md).

## What this does *not* do (yet)

- **No vector-field functionals.** `evaluate(Field<Vec3d>, …)` arrives in
  Phase 11+ when the higher-rank functionals first need it.
- **No higher-order derivatives in `f`.** Phase 6 adds biharmonic and
  mixed partials; Phase 11 adds up-to-4th-order gradients.
- **No variational derivative `δF/δψ`.** Phase 12 (CH demo) adds
  closed-form variational derivatives where they are needed.
- **No fused-loop optimization.** `evaluate` materializes intermediate
  Fields for `∇ψ` and/or `∇²ψ`. Phase 11+ replaces this with
  expression-template fusion.

## Quick reference: which arity to pick

- Integrand depends only on `ψ` (e.g., `∫ ψ² dV`, `∫ W(ψ²-1)² dV`):
  use `f(double)`. Saves both stencil computations.
- Integrand uses `∇ψ` (e.g., the gradient term of any GL-style energy):
  use `f(double, const Vec3d&)`. Saves the Laplacian.
- Integrand uses `∇²ψ` (Cahn–Hilliard kernel, surface-tension
  derivatives): use `f(double, const Vec3d&, double)`.
