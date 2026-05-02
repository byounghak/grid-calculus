# User Guide — Phase 3: Domain integration

> Spec dir: `specs/2026-05-02-phase-3-domain-integration/`. Version: `0.3.0`.

## What this gives you

A primitive that turns a scalar `Field<double>` into a single number — the
discrete integral

$$
\int_{\Omega} f(\mathbf{x})\, dV \;\approx\; h_x h_y h_z \sum_{i,j,k} f_{ijk}
$$

over the periodic Cartesian domain `[0, Nx*hx) × [0, Ny*hy) × [0, Nz*hz)`.

This is the integration backbone the upcoming functional-evaluation API will
sit on (Phase 4).

## API

Header: [`include/gridcalc/func/Integrate.hpp`](../../../include/gridcalc/func/Integrate.hpp).

```cpp
namespace gridcalc::func {

struct Pairwise {};   // tag selecting pairwise recursive summation
struct Kahan {};      // tag selecting Kahan/Neumaier compensated summation

double integrate(const core::Field<double>& field, Pairwise = {}) noexcept;
double integrate(const core::Field<double>& field, Kahan)        noexcept;

}  // namespace gridcalc::func
```

Pick the reducer with a tag:

- `integrate(f)` or `integrate(f, Pairwise{})` — **default**. Pairwise
  recursive summation; rounding error grows like `O(ε log n)`. Cheap. Use
  this unless you have a reason not to.
- `integrate(f, Kahan{})` — Neumaier-style compensated summation;
  ~bit-perfect accuracy at ~3-4× cost per element. Use when integrating
  long sums of values whose magnitudes vary by many orders of magnitude.

The grid's cell volume `hx * hy * hz` is read from the input `Field`'s
`Grid` and applied by `integrate()` itself. You do not multiply by `dV`
yourself.

Periodic wrap is irrelevant here: every grid point is sampled exactly
once, regardless of the input field's `IndexPolicy`.

## Worked example

Volume of a 2π cube:

```cpp
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/func/Integrate.hpp>

const int N = 32;
const double L = 2.0 * M_PI;
const double h = L / N;

gridcalc::core::Grid grid(N, N, N, gridcalc::core::Vec3d(h, h, h));
gridcalc::core::Field<double> ones(grid, 1.0);

double V = gridcalc::func::integrate(ones);
// V == L*L*L == (2*pi)^3, to machine precision.
```

A periodic sinusoid integrates to zero (every axis-orthogonal trig
component sums exactly to zero on equispaced periodic samples):

```cpp
gridcalc::core::Field<double> psi(grid, [](double x, double y, double z) {
  return std::cos(x) * std::sin(2.0 * y) * std::sin(3.0 * z);
});
double I = gridcalc::func::integrate(psi);   // I ~ 0  (within ~1e-15)
```

## What this does *not* do (yet)

- **No callable-integrand overload.** `∫ f(ψ, ∇ψ, ∇²ψ) dV` for
  user-supplied `f` is the headline of Phase 4 (functional evaluation).
- **No vector-field overload.** `integrate(Field<Vec3d>)` arrives when the
  Phase 11+ higher-rank functionals first need it.
- **No threading.** Phase 3 is single-threaded; Phase 20 (performance
  pass) parallelizes the reducer while preserving bit-identical output.
- **No non-periodic boundary handling.** The midpoint formula assumes the
  domain is periodic. Once Phase 18 introduces Dirichlet/Neumann fields,
  endpoint corrections (Euler–Maclaurin, Simpson) become relevant.

For the why-this-formula-is-correct background — derivation, error
analysis, references — see
[`docs/developer-note/notes/phase-3-domain-integration.md`](../../developer-note/notes/phase-3-domain-integration.md).
