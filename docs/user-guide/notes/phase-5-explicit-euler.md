# User Guide — Phase 5: Explicit Euler diffusion solver

> Spec dir: `specs/2026-05-02-phase-5-explicit-euler/`. Version: `0.5.0`.

## What this gives you

Time integration enters the picture. Two new public functions in the
`gridcalc::solve` namespace:

- **`explicitEuler`** — a generic forward-Euler integrator over an
  arbitrary RHS callable.
- **`diffuse`** — a thin convenience driver for the heat equation
  $\partial_t \psi = D \nabla^2 \psi$, with a built-in CFL stability check.

Both mutate the input `Field<double>` **in place** — the natural pattern
for iterative time stepping. If you want to keep the initial condition,
copy the field first.

## API

Headers:
[`include/gridcalc/solve/ExplicitEuler.hpp`](../../../include/gridcalc/solve/ExplicitEuler.hpp),
[`include/gridcalc/solve/Diffusion.hpp`](../../../include/gridcalc/solve/Diffusion.hpp).

```cpp
namespace gridcalc::solve {

template <class Rhs>
void explicitEuler(core::Field<double>& psi,
                   Rhs&& rhs,
                   double dt,
                   int n_steps);

void diffuse(core::Field<double>& psi,
             double D,
             double dt,
             int n_steps);

}  // namespace gridcalc::solve
```

### `explicitEuler` semantics

At each step: `tmp = rhs(psi); psi[n] += dt * tmp[n]`. The RHS callable
must be invocable as `Field<double> rhs(const Field<double>&)`; mismatch
trips a clean `static_assert` at the call site.

The integrator does **not** validate stability — that depends on the
specific RHS, which is invisible to the integrator. If your RHS has a
known stability bound, check it yourself before calling.

Throws `std::invalid_argument` on `dt < 0` or `n_steps < 0`.

### `diffuse` semantics

Advances the heat equation by `n_steps` explicit-Euler steps of size
`dt`. Internally calls `explicitEuler` with `diff::laplacian` as the
RHS and `D * dt` as the effective step (mathematically identical to
running an unscaled Laplacian RHS at a scaled time step; the trick
avoids a per-step scaling loop).

**CFL stability check.** Before integration begins, `diffuse` validates
the von Neumann bound

$$
D \cdot dt \cdot \left( \frac{1}{h_x^2} + \frac{1}{h_y^2} + \frac{1}{h_z^2} \right) \;\le\; \frac{1}{2}
$$

and throws `std::invalid_argument` (with a message naming the offending
values) if the bound is exceeded. For isotropic spacing this reduces to
$D\,dt/h^2 \le 1/(2d) = 1/6$ in 3D, the form quoted in the roadmap.

Throws on negative `D`, `dt`, or `n_steps`.

## Worked example

Diffuse a trig eigenfunction and compare to the analytical decay:

```cpp
#include <cmath>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/solve/Diffusion.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
using gridcalc::solve::diffuse;

constexpr double kPi = 3.14159265358979323846;

const int N = 32;
const double L = 2.0 * kPi;
const double h = L / N;
const double D = 0.1;
const double dt = 0.05;
const int n_steps = 100;
const double T = n_steps * dt;  // 5.0

Grid grid(N, N, N, Vec3d(h, h, h));

// Initial condition: psi_0 = sin(x) sin(y) sin(z)
// is an eigenfunction of the Laplacian with eigenvalue -3.
Field<double> psi(grid, [](double x, double y, double z) {
  return std::sin(x) * std::sin(y) * std::sin(z);
});

// In place: psi is mutated.
diffuse(psi, D, dt, n_steps);

// After integration, psi ~= exp(-3 D T) * sin(x) sin(y) sin(z),
// to within stencil-order error (~few percent at N=32).
```

## Custom RHS via `explicitEuler`

If you have an RHS that isn't $D\nabla^2\psi$, call `explicitEuler`
directly and provide your own callable. **You** are responsible for
the relevant stability check.

```cpp
#include <gridcalc/diff/Laplacian.hpp>
#include <gridcalc/solve/ExplicitEuler.hpp>

// Reaction-diffusion: d_t psi = D laplacian(psi) - k * psi
using gridcalc::diff::laplacian;
using gridcalc::solve::explicitEuler;

const double D = 0.1;
const double k = 0.05;

auto rhs = [D, k](const Field<double>& f) {
  Field<double> out = laplacian(f);
  // out *= D; subtract k * f
  const std::size_t N = out.getSize();
  double* o = out.data();
  const double* p = f.data();
  for (std::size_t n = 0; n < N; ++n) {
    o[n] = D * o[n] - k * p[n];
  }
  return out;
};

explicitEuler(psi, rhs, dt, n_steps);
```

## What this does *not* do (yet)

- **No higher-order time integrators.** RK4 + a generic `solve::Integrator`
  interface arrive in Phase 6; the explicit-Euler entry point will be
  refactored to consume that interface then.
- **No adaptive time stepping.**
- **No implicit / Crank–Nicolson diffusion.** Phase 19 introduces
  implicit diffusion, after Phase 18 lands non-periodic boundary policies.
- **No automatic CFL adjustment.** `diffuse` throws; the caller picks a
  smaller `dt` and retries.
- **No threading.** Phase 20 (performance pass) parallelizes the
  per-step accumulation loop.

For the why-this-stability-bound and the time-discretization error
analysis, see
[`docs/developer-note/notes/phase-5-explicit-euler.md`](../../developer-note/notes/phase-5-explicit-euler.md).
