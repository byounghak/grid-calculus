# User Guide — Phase 6: RK4 + generic time integrator interface

> Spec dir: `specs/2026-05-02-phase-6-rk4-integrator/`. Version: `0.6.0`.

## What this gives you

A higher-order time integrator (classic 4-stage Runge–Kutta) and a
**generic dispatcher** that lets you pick which integrator to use:

```cpp
solve::integrate(psi, rhs, dt, n_steps, IntegratorTag{});
```

with `IntegratorTag` being either `solve::ExplicitEuler` (Phase 5) or
`solve::RK4` (new in Phase 6).

`solve::diffuse` is now templated on the integrator tag so you can swap
between Euler and RK4 with one extra argument:

```cpp
solve::diffuse(psi, D, dt, n_steps);                    // Euler (default)
solve::diffuse(psi, D, dt, n_steps, solve::RK4{});      // RK4
```

RK4 is **more accurate** at the same `dt` (4 orders better — the
time-error per step is `O(dt^5)` vs Euler's `O(dt^2)`) and has a
**larger stability bound** for the heat equation (`D dt sum(1/h_a^2) ≤ 0.6963`
vs Euler's `0.5`). The cost is 4 RHS evaluations per step instead of 1.

## API

Headers:
[`include/gridcalc/solve/Integrator.hpp`](../../../include/gridcalc/solve/Integrator.hpp)
(one-stop aggregator),
[`include/gridcalc/solve/ExplicitEuler.hpp`](../../../include/gridcalc/solve/ExplicitEuler.hpp),
[`include/gridcalc/solve/RK4.hpp`](../../../include/gridcalc/solve/RK4.hpp),
[`include/gridcalc/solve/Diffusion.hpp`](../../../include/gridcalc/solve/Diffusion.hpp).

```cpp
namespace gridcalc::solve {

// Integrator tag types -- empty structs that carry compile-time properties.
struct ExplicitEuler { static constexpr double diffusionCFLLimit = 0.5;    };
struct RK4           { static constexpr double diffusionCFLLimit = 0.6963; };

// Generic dispatcher (overloaded on tag).
template <class Rhs>
void integrate(core::Field<double>& psi, Rhs&& rhs, double dt, int n_steps, ExplicitEuler);

template <class Rhs>
void integrate(core::Field<double>& psi, Rhs&& rhs, double dt, int n_steps, RK4);

// Phase-5 thin wrapper -- equivalent to integrate(..., ExplicitEuler{}).
template <class Rhs>
void explicitEuler(core::Field<double>& psi, Rhs&& rhs, double dt, int n_steps);

// Diffusion driver -- now templated on integrator tag.
template <class Integrator = ExplicitEuler>
void diffuse(core::Field<double>& psi, double D, double dt, int n_steps,
             Integrator tag = {});

}  // namespace gridcalc::solve
```

### Picking an integrator

| Want                                                  | Use         |
|-------------------------------------------------------|-------------|
| Simplest, smallest per-step cost                      | `ExplicitEuler{}` (Phase 5 default) |
| Best accuracy per step                                | `RK4{}` |
| Larger time step (closer to RK4's `~1.39×` looser CFL) | `RK4{}` |
| Pure linear / Hamiltonian / oscillatory dynamics      | `RK4{}` (lower phase error) |
| One-pass throwaway integration of a smooth IC         | Either; default Euler is fine |

Both integrators take the same RHS callable signature:
`Field<double> rhs(const Field<double>&)`. SFINAE-validates; mismatch
trips a `static_assert`.

### CFL stability check

`solve::diffuse` validates the von-Neumann bound
`D · dt · sum_a(1/h_a^2) ≤ Tag::diffusionCFLLimit` before integration
begins and throws `std::invalid_argument` if violated. The error
message names the integrator type, the offending values, and the
applicable limit.

`solve::integrate` and `solve::explicitEuler` do **not** check CFL —
they're generic over the RHS, which means they cannot know the
relevant stability bound. If you call them directly with a custom
RHS, you own the stability argument.

## Worked example: heat equation with RK4

```cpp
#include <cmath>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/solve/Diffusion.hpp>
#include <gridcalc/solve/Integrator.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
using gridcalc::solve::diffuse;
using gridcalc::solve::RK4;

constexpr double kPi = 3.14159265358979323846;

const int N = 32;
const double L = 2.0 * kPi;
const double h = L / N;
const double D = 0.1;
const double dt = 0.05;
const int n_steps = 100;

Grid grid(N, N, N, Vec3d(h, h, h));

Field<double> psi(grid, [](double x, double y, double z) {
  return std::sin(x) * std::sin(y) * std::sin(z);
});

// One-line switch from Euler to RK4: pass the integrator tag.
diffuse(psi, D, dt, n_steps, RK4{});

// psi now matches the analytical exp(-3 D T) sin(x) sin(y) sin(z) to
// within ~0.5% relative max-norm error -- 10x tighter than Phase 5's
// Euler at the same parameters.
```

## Worked example: custom RHS (no CFL gating)

```cpp
#include <gridcalc/diff/Laplacian.hpp>
#include <gridcalc/solve/Integrator.hpp>

// Reaction-diffusion: dpsi/dt = D laplacian(psi) - k psi
auto rhs = [D, k](const Field<double>& f) {
  Field<double> out = gridcalc::diff::laplacian(f);
  const std::size_t N_ = out.getSize();
  double* o = out.data();
  const double* p = f.data();
  for (std::size_t n = 0; n < N_; ++n) {
    o[n] = D * o[n] - k * p[n];
  }
  return out;
};

// RK4 makes this much more accurate than Euler at the same dt.
gridcalc::solve::integrate(psi, rhs, dt, n_steps, gridcalc::solve::RK4{});
```

## What this does *not* do (yet)

- **No adaptive time stepping.** Embedded RK methods (RK4(5), Dormand–Prince)
  arrive only if a downstream phase needs them.
- **No symplectic / structure-preserving integrators.** Heat equation is
  dissipative; symplectic schemes are for Hamiltonian flows.
- **No multistep methods (Adams–Bashforth, BDF).** Implicit BDF arrives in
  Phase 19.
- **No persistent scratch pool for RK4.** Each step allocates ~6 fields;
  Phase 20 swaps in a reusable scratch buffer.
- **Boundary conditions: still periodic-only.** Phase 18's Dirichlet/Neumann
  policies will need their own integrator-test coverage; tracked in
  `STATUS.md` "Open / deferred items."

For the why-RK4-is-order-4 derivation and the stability-bound calculation,
see [`docs/developer-note/notes/phase-6-rk4-integrator.md`](../../developer-note/notes/phase-6-rk4-integrator.md).
