# User Guide — Phase 1: Periodic 3D scalar grid + Laplacian

> Spec dir: `specs/2026-05-01-phase-1-laplacian/`. Version: `0.1.0`.
> *Note: this note was written retroactively after Phase 4, so the
> "current PR" framing of later notes is absent here.*

## What this gives you

The structural foundation of the library: a Cartesian-orthogonal **sampling
grid**, a **periodic scalar field** built on it, and a **2nd-order
Laplacian** operator. Everything later in the roadmap (gradient,
divergence, integration, functional evaluation, …) sits on top of these
three pieces.

## API

Headers:
[`include/gridcalc/core/Grid.hpp`](../../../include/gridcalc/core/Grid.hpp),
[`include/gridcalc/core/Field.hpp`](../../../include/gridcalc/core/Field.hpp),
[`include/gridcalc/core/IndexPolicy.hpp`](../../../include/gridcalc/core/IndexPolicy.hpp),
[`include/gridcalc/stencil/CentralDifference.hpp`](../../../include/gridcalc/stencil/CentralDifference.hpp),
[`include/gridcalc/diff/Laplacian.hpp`](../../../include/gridcalc/diff/Laplacian.hpp).

```cpp
namespace gridcalc::core {

class Grid {
 public:
  Grid(int Nx, int Ny, int Nz, const Vec3d& cell_size);
  int          getNx() const noexcept;
  int          getNy() const noexcept;
  int          getNz() const noexcept;
  const Vec3d& getCellSize() const noexcept;
  double       getCellVolume() const noexcept;
  std::size_t  getSize() const noexcept;
  std::size_t  getLinearIndex(int i, int j, int k) const noexcept;
};

template <class T, class Policy = IndexPolicy::Periodic>
class Field {
 public:
  explicit Field(const Grid& grid);                // zero-initialized
  Field(const Grid& grid, const T& value);         // broadcast
  template <class F> Field(const Grid& grid, F&& f);  // sample f(x, y, z)

  T&       operator()(int i, int j, int k) noexcept;     // wraps via Policy
  const T& operator()(int i, int j, int k) const noexcept;

  const Grid&  getGrid() const noexcept;
  std::size_t  getSize() const noexcept;
  T*           data() noexcept;
  const T*     data() const noexcept;
  T*           begin() noexcept;  T* end() noexcept;
};

namespace IndexPolicy {
  struct Periodic { static constexpr int wrap(int i, int N) noexcept; };
}

}  // namespace gridcalc::core

namespace gridcalc::diff {
core::Field<double> laplacian(const core::Field<double>& field);
}
```

### Storage layout

`Field<T>` storage is **i-fastest**: the linear offset for `(i, j, k)` is
`i + Nx * (j + Ny * k)`. This matches the Eigen / FFT convention and is
part of the public contract from Phase 1 forward — Phase 9 (spectral
verification) and Phase 20 (benchmarks) both rely on it.

### Index policy

`Field<T, Policy>` defaults to `IndexPolicy::Periodic`, which wraps
out-of-range indices modulo `N`. `field(-1, j, k)` returns
`field(Nx - 1, j, k)`. `Dirichlet` and `Neumann` policies arrive in
Phase 18; the policy template parameter is the only piece of
forward-looking design admitted at Phase 1.

## Worked example

Recover the Laplacian eigenvalue of a trig product:

```cpp
#include <cmath>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/diff/Laplacian.hpp>

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
using gridcalc::diff::laplacian;

constexpr double kPi = 3.14159265358979323846;

const int N = 32;
const double L = 2.0 * kPi;
const double h = L / N;

Grid grid(N, N, N, Vec3d(h, h, h));
Field<double> psi(grid, [](double x, double y, double z) {
  return std::sin(x) * std::sin(y) * std::sin(z);
});

Field<double> lap = laplacian(psi);
// lap ~= -3 * psi  (since the analytical Laplacian eigenvalue is -3)
// Relative max-norm error ~ k^4 h^2 / (12 k^2) ~ h^2 / 12 ~ 3e-3 at N = 32.
```

## What this does *not* do (yet)

- **Periodic boundary only.** Dirichlet / Neumann policies arrive in Phase 18.
- **2nd-order accuracy only.** Higher-order accuracy (`Order = 4`) on the
  same operator arrives in Phase 7.
- **Scalar `T = double` only in practice.** Although `Field<T>` is
  templated, only `double` is exercised at Phase 1. `Field<Vec3d>` arrives
  in Phase 2's gradient output.

For the why-this-stencil-converges background — derivation, error
analysis, references — see
[`docs/developer-note/notes/phase-1-laplacian.md`](../../developer-note/notes/phase-1-laplacian.md).
