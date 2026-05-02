# Plan: Phase 1 — Periodic 3D scalar grid + Laplacian

Five commits on branch `2026-05-01-phase-1-laplacian`. Order keeps the build green between commits; specs land first, then the foundation types in dependency order, then the operator, then the version bump.

## Group 1 — Phase 1 spec files

- Add `specs/2026-05-01-phase-1-laplacian/{plan,requirements,validation}.md`.

**Commit message:** `docs: add Phase 1 specs (plan, requirements, validation)`

## Group 2 — IndexPolicy + Grid

- `include/gridcalc/core/IndexPolicy.hpp`:
  - Namespace `gridcalc::core::IndexPolicy`.
  - `struct Periodic` with `static constexpr int wrap(int i, int N) noexcept`. Document `\since 0.1.0`.
- `include/gridcalc/core/Grid.hpp`:
  - Namespace `gridcalc::core`.
  - `class Grid` carrying `int _Nx, _Ny, _Nz` and `Vec3d _cell_size` (where `Vec3d = Eigen::Vector3d`).
  - Constructor: `Grid(int Nx, int Ny, int Nz, const Vec3d& cell_size)`. Throws `std::invalid_argument` on `Nx <= 0`, `Ny <= 0`, `Nz <= 0`, or any non-positive `cell_size` component.
  - Accessors: `getNx`, `getNy`, `getNz`, `getCellSize`, `getCellVolume`, `getSize`, `size` (STL shim → `getSize`), `getLinearIndex(int i, int j, int k) -> std::size_t` (assumes valid in-range; computes `i + Nx*(j + Ny*k)`).
  - All methods inline; header-only.
  - Pull `Vec3d` from a small new `include/gridcalc/core/EigenAliases.hpp` (or define inline here — pick the inline option to keep Phase 1 minimal; promote to a shared header at Phase 2 when `Vec3d` is needed by Gradient).
  - Doxygen on every public member; `\since 0.1.0`.
- `test/grid_test.cpp`:
  - `GridTest.ConstructsAndExposesShape` — values match constructor inputs.
  - `GridTest.ComputesCellVolume` — `hx*hy*hz`.
  - `GridTest.LinearIndexFormula` — spot-check a few `(i,j,k)` triples against the i-fastest formula.
  - `GridTest.RejectsZeroOrNegativeShape` and `GridTest.RejectsNonPositiveSpacing` — death/exception tests.
- `test/index_policy_test.cpp`:
  - `PeriodicTest.WrapsInRange` — `wrap(0..N-1, N) == self`.
  - `PeriodicTest.WrapsNegative` — `wrap(-1, 5) == 4`, `wrap(-7, 5) == 3`.
  - `PeriodicTest.WrapsLarge` — `wrap(5, 5) == 0`, `wrap(13, 5) == 3`.
- `test/CMakeLists.txt` — add the two new sources; reuse existing `gridcalc_tests` executable.

**Commit message:** `feat(core): add Grid and Periodic IndexPolicy`

## Group 3 — Field

- `include/gridcalc/core/Field.hpp`:
  - `template <class T, class Policy = IndexPolicy::Periodic> class Field`.
  - Members: `Grid _grid`, `std::vector<T> _data`.
  - Constructors: zero-init, value-init, callable-init (callable signature `T(double, double, double)`).
  - `T& operator()(int i, int j, int k)` and const overload — apply `Policy::wrap` per axis, then `_grid.getLinearIndex(...)` for storage offset.
  - Accessors: `getGrid`, `getSize`, `size` (STL shim), `data`, `begin`, `end` (and const variants).
  - Inline; header-only.
  - Doxygen + `\since 0.1.0`.
- `test/field_test.cpp`:
  - `FieldTest.ZeroInitFillsZero`.
  - `FieldTest.ValueInitBroadcasts`.
  - `FieldTest.CallableInitSamplesCartesian` — verify a known function (e.g. `(x,y,z) -> x + 2y + 3z`) at a few points.
  - `FieldTest.PeriodicWrapAtNegativeIndex` — `field(-1, 0, 0) == field(Nx-1, 0, 0)`.
  - `FieldTest.PeriodicWrapAtOverflowIndex` — `field(Nx, 0, 0) == field(0, 0, 0)`.
  - `FieldTest.WriteThenReadRoundTrip`.
- `test/CMakeLists.txt` — add the new source.

**Commit message:** `feat(core): add Field<T, Policy=Periodic> with three constructors`

## Group 4 — CentralDifference + Laplacian + acceptance tests

- `include/gridcalc/stencil/CentralDifference.hpp`:
  - Namespace `gridcalc::stencil`.
  - `template <int Order> struct Coefficients;` (primary, undefined).
  - `template <> struct Coefficients<2>` exposing `static constexpr int radius = 1`, `static constexpr std::array<int, 3> offsets = {-1, 0, 1}`, `static constexpr std::array<double, 3> weights = {1.0, -2.0, 1.0}`. Document that weights are *unscaled* — the consumer multiplies by `1/h^2`.
- `include/gridcalc/diff/Laplacian.hpp`:
  - Namespace `gridcalc::diff`.
  - `core::Field<double> laplacian(const core::Field<double>& field)`.
  - Implementation: triple loop over `(i, j, k) in [0, Nx) x [0, Ny) x [0, Nz)`. For each axis, sum `weights[m] * field(i + offsets[m], j, k)` (etc.) and divide by `h_axis^2`. Sum the three axis contributions.
  - Inline.
- `test/laplacian_test.cpp`:
  - `LaplacianTest.RecoversTrigEigenvalue` — $L = 2\pi$, $N = 32$, $\mathbf{k} = (1, 2, 1)$. Build `Field` from $\sin(k_x x)\sin(k_y y)\sin(k_z z)$ via the callable constructor; compute `laplacian(f)`; verify `laplacian / f ≈ -(kx^2 + ky^2 + kz^2) = -6` to relative max-norm tolerance `5e-3`.
  - `LaplacianTest.SecondOrderConvergence` — sweep $N \in \{16, 32, 64\}$ on $L = 2\pi$, $\mathbf{k} = (1, 1, 1)$. Compute relative max-norm error at each resolution. Assert error decreases by factor in `[3.5, 4.5]` per doubling; assert least-squares log-log slope in `[1.8, 2.2]`.
  - `LaplacianTest.OutputGridMatchesInput` — sanity: `laplacian(f).getGrid().getSize() == f.getGrid().getSize()`.
- `test/CMakeLists.txt` — add the new source.

**Commit message:** `feat(diff): add 2nd-order periodic Laplacian + convergence test`

## Group 5 — Version bump

- `CMakeLists.txt`: bump `project(gridcalc VERSION 0.0.1 ...)` → `0.1.0`. Per `tech-stack.md`'s pre-1.0 policy, minor bumps may break API; a real public surface (Grid, Field, Laplacian) now exists, so 0.1.0 is the right level.
- `test/version_test.cpp`: update the expected string from `"0.0.1"` to `"0.1.0"`.
- `specs/STATUS.md`: update **Where the project stands**, the **Next action** block (point at Phase 2), and bump the Phase 1 row to "Done" — landed as the same commit so STATUS doesn't lag the merge.

**Commit message:** `build: bump version to 0.1.0 and mark Phase 1 done`
