# User Guide — Phase 2: Gradient and divergence (scalar)

> Spec dir: `specs/2026-05-02-phase-2-gradient-divergence/`. Version: `0.2.0`.
> *Note: this note was written retroactively after Phase 4 (see
> specs/STATUS.md "Decisions worth knowing").*

## What this gives you

The basic first-order vector operators on a periodic scalar field:

- **Gradient:** `gradient(ψ) -> Field<Vec3d>`. At each grid point, returns
  the Cartesian triple $(\partial_x \psi,\, \partial_y \psi,\, \partial_z \psi)$.
- **Divergence:** `divergence(V) -> Field<double>`. For a vector field
  $\mathbf{V}$ stored as `Field<Vec3d>`, returns
  $\partial_x V_x + \partial_y V_y + \partial_z V_z$ at each grid point.

Both are 2nd-order accurate central differences applied per axis, with
periodic wrap delegated to the input `Field`'s `Policy` (`Periodic` from
Phase 1 onward).

## API

Headers:
[`include/gridcalc/diff/Gradient.hpp`](../../../include/gridcalc/diff/Gradient.hpp),
[`include/gridcalc/diff/Divergence.hpp`](../../../include/gridcalc/diff/Divergence.hpp),
[`include/gridcalc/stencil/FirstDerivative.hpp`](../../../include/gridcalc/stencil/FirstDerivative.hpp),
[`include/gridcalc/core/EigenAliases.hpp`](../../../include/gridcalc/core/EigenAliases.hpp).

```cpp
namespace gridcalc::core {
using Vec3d = Eigen::Vector3d;
}

namespace gridcalc::diff {

core::Field<core::Vec3d> gradient   (const core::Field<double>&)         noexcept;
core::Field<double>      divergence (const core::Field<core::Vec3d>&)    noexcept;

}  // namespace gridcalc::diff
```

The `Vec3d` alias was relocated from `core/Grid.hpp` (Phase 1) to
`core/EigenAliases.hpp` in Phase 2 once a second operator needed it.
The fully-qualified name `gridcalc::core::Vec3d` is unchanged from
Phase 1.

## Worked example

Compute a gradient, then a divergence, then verify that
$\nabla \cdot \nabla \psi$ converges to $\nabla^2 \psi$ at order 2:

```cpp
#include <cmath>
#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/diff/Divergence.hpp>
#include <gridcalc/diff/Gradient.hpp>
#include <gridcalc/diff/Laplacian.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
using gridcalc::diff::divergence;
using gridcalc::diff::gradient;
using gridcalc::diff::laplacian;

constexpr double kPi = 3.14159265358979323846;

const int N = 32;
const double L = 2.0 * kPi;
const double h = L / N;

Grid grid(N, N, N, Vec3d(h, h, h));
Field<double> psi(grid, [](double x, double y, double z) {
  return std::sin(x) * std::sin(y) * std::sin(z);
});

Field<Vec3d>  grad_psi  = gradient(psi);
Field<double> lap_psi   = laplacian(psi);
Field<double> indirect  = divergence(grad_psi);
// indirect ~ lap_psi at order 2 in h, but the two are *not* identically
// equal at finite h -- divergence(gradient(psi)) is a stride-2 stencil
// while laplacian(psi) is a 3-point stencil. They share the same
// continuum limit but have different leading-order error constants.
```

## Vector-field example for divergence

`divergence` takes any `Field<Vec3d>`, not only one produced by
`gradient`. For $\mathbf{V} = (\sin x \cos y,\, \cos x \sin y,\, \sin z)$
on $[0, 2\pi]^3$, the analytical divergence is
$2\cos x \cos y + \cos z$:

```cpp
Field<Vec3d> V(grid, [](double x, double y, double z) {
  return Vec3d(std::sin(x) * std::cos(y),
               std::cos(x) * std::sin(y),
               std::sin(z));
});
Field<double> div_V = divergence(V);
// div_V ~ 2 cos(x) cos(y) + cos(z), with relative max-norm error O(h^2).
```

## What this does *not* do (yet)

- **No curl operator.** Vector calculus on `Field<Vec3d>` beyond
  divergence arrives in later phases.
- **No 4th-order accuracy.** `stencil::FirstDerivative<4>` is deferred
  to Phase 7 alongside `Coefficients<4>`.
- **Only periodic wrap.** Dirichlet / Neumann gradients arrive in
  Phase 18.
- **No per-component-decomposition helper.** `Field<Vec3d>` is the
  canonical gradient return; iterate over it and project per axis if you
  need a `Field<double>` for one component.

For the why-divergence(gradient)-isn't-laplacian background and the
full error analysis — see
[`docs/developer-note/notes/phase-2-gradient-divergence.md`](../../developer-note/notes/phase-2-gradient-divergence.md).
