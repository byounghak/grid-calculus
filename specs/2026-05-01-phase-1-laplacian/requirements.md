# Requirements: Phase 1 — Periodic 3D scalar grid + Laplacian

## Roadmap reference

> **Goal.** First real numerics: a scalar field on a periodic 3D grid with a 2nd-order-accurate Laplacian.
>
> **Deliverables.**
>
> - `core/Grid.hpp` — Nx/Ny/Nz, lattice vectors, cell volume. Cartesian-orthogonal only for now.
> - `core/Field.hpp` — `Field<double>` with linear indexing.
> - `core/IndexPolicy.hpp` — `Periodic` policy with wrapping `index(i,j,k)`.
> - `stencil/CentralDifference.hpp` — 2nd-order central differences in 1D.
> - `diff/Laplacian.hpp` — `laplacian(field) -> Field<double>`.
>
> **Acceptance.**
>
> - Apply Laplacian to $\sin(k_x x)\sin(k_y y)\sin(k_z z)$ → recover $-(k_x^2+k_y^2+k_z^2)$ times the input.
> - Convergence test: halving $h$ reduces error by 4×, log-log slope ≈ 2.

Refer to [mission.md](../mission.md), [tech-stack.md](../tech-stack.md) (especially the code-style section), and [roadmap.md](../roadmap.md) for global context. The decision recorded in STATUS.md to defer `byounghak/atomic-structure` integration to Phase 15+ applies here: Phase 1 `Grid` is an in-house Cartesian struct, not a wrapper over `NSPC_Atomic_Structure::Structure`.

## Decisions made for this feature

- **Grid cell size: per-axis `Vec3d {hx, hy, hz}`** stored as an `Eigen::Vector3d` member. Forward-compatible to anisotropic Cartesian grids without an API break. Constructor signature: `Grid(int Nx, int Ny, int Nz, const Vec3d& cell_size)`.
- **`Field<T>` storage layout: i fastest** — `linear_index = i + Nx * (j + Ny * k)`. Matches FFT libraries (PocketFFT, FFTW) and Eigen's column-major default; cache-friendly for x-direction stencil sweeps. The layout is part of the public contract from Phase 1 forward (benchmarks at Phase 20 will reference it).
- **`Field<T>` construction overloads (three):**
  1. `Field(const Grid& grid)` — zero-initialized.
  2. `Field(const Grid& grid, T value)` — broadcast a single value.
  3. `template <class F> Field(const Grid& grid, F&& f)` — Cartesian sample of a callable `f(double x, double y, double z) -> T`. Sample point `(i, j, k)` evaluates at `(i*hx, j*hy, k*hz)` (origin at the grid corner; no half-cell offset).
- **`Field<T>` storage policy: by-value `Grid`.** `Field` holds `Grid _grid` by value (Grid is ~36 bytes; a copy is cheaper than a lifetime-management contract). Two `Field`s built from the "same" `Grid` carry equal-but-distinct copies; `Field::getGrid()` returns const ref to the internal copy.
- **`IndexPolicy::Periodic`** is a stateless struct with a `static constexpr int wrap(int i, int N) noexcept` returning `((i % N) + N) % N`. `Field<T>::operator()(int i, int j, int k)` applies the wrap on every access (Phase 1 prioritizes ergonomics; the branch is predictable and Phase 20 can revisit if it shows up in profiles).
- **`Field<T>` is templated on Policy with default `Periodic`** — `template <class T, class Policy = IndexPolicy::Periodic> class Field`. Pre-emptively makes Phase 18 (Dirichlet/Neumann) a non-breaking addition; the cost at Phase 1 is one template parameter and one default. This is the only forward-looking template parameter we're introducing — everything else stays YAGNI.
- **`stencil/CentralDifference.hpp`: `Coefficients<Order>` template, `Order=2` specialization only.** Phase 7 adds the `Order=4` specialization with no refactor of call sites or `Laplacian`. The `Order=2` specialization exposes `static constexpr` arrays of integer offsets and double weights for the second derivative ($f_{i-1} - 2 f_i + f_{i+1}$, *unscaled* — the `1/h^2` factor is applied by the consumer). Stencil radius is exposed as `static constexpr int radius`.
- **`diff::laplacian` signature: allocating only.** `core::Field<double> laplacian(const core::Field<double>& field)` matches the roadmap text. An in-place / out-parameter overload is deferred to Phase 5 when time-stepping needs it (avoids a Phase-1 API decision that Phase 5 may revise).
- **Free-function naming.** Snake_case for non-method free functions (`laplacian`, future `gradient`, `divergence`). Matches the roadmap deliverable text. The code-style policy in `tech-stack.md` only pins method naming; this adds the free-function convention by precedent for future phases.
- **Acceptance-test inputs.**
  - **Eigenvalue test:** $L = 2\pi$, $N_x = N_y = N_z = 32$, $\mathbf{k} = (1, 2, 1)$. Eigenvalue $-(k_x^2 + k_y^2 + k_z^2) = -6$. Tolerance: relative max-norm error $< 1.5\times 10^{-2}$. The leading-order error of the 2nd-order central-difference eigenvalue on a trig product is $(k_x^4+k_y^4+k_z^4)\,h^2 / [12 (k_x^2+k_y^2+k_z^2)]$; for $\mathbf{k}=(1,2,1)$ at $h=2\pi/32$ that evaluates to $\approx 9.6\times 10^{-3}$, so $1.5\times 10^{-2}$ leaves ~50% headroom.
  - **Convergence test:** $L = 2\pi$, $\mathbf{k} = (1, 1, 1)$, sweep $N \in \{16, 32, 64\}$. Compute relative max-norm error at each $N$; assert error decreases by a factor in $[3.5, 4.5]$ between successive resolutions; assert least-squares log-log slope of error vs $h$ in $[1.8, 2.2]$.
- **Header-only at Phase 1.** All Phase-1 code lives in headers under `include/gridcalc/`. The `gridcalc` CMake target stays `INTERFACE`. Promoting to a real STATIC library is deferred to the first phase that has compilation work that doesn't belong inline (likely Phase 5 or 9).

## Non-goals

- **No half-cell-offset / staggered grids.** Sample point $(i, j, k)$ sits exactly at $(i\,h_x, j\,h_y, k\,h_z)$ from the origin.
- **No non-uniform grids.** Spacing is constant per axis.
- **No higher-order accuracy** at Phase 1. 4th-order central differences land in Phase 7.
- **No gradient or divergence operators.** Both land in Phase 2.
- **No FFT comparison.** That is Phase 9's CI gate.
- **No OpenMP.** Single-threaded baseline; threading lands in Phase 20 alongside the benchmark suite.
- **No in-place / out-parameter `laplacian` overload.** Deferred to Phase 5.
- **No `Field<Vector3d>` / non-scalar fields.** Land in Phase 13.
- **No `\since` lint enforcement** — that is Phase 10's job. Phase 1 carries `\since 0.1.0` tags on every new public symbol so the future lint passes when it lands.

## Deferred questions

- **Promotion of `gridcalc` from INTERFACE to STATIC** — when the first non-inline compilation unit shows up.
- **`Field<T>::operator()` wrap-on-every-access vs. caller-side wrap** — Phase 20 (performance pass) decides whether to keep the unconditional wrap or introduce an unchecked accessor.
- **Pre-release tag for Phase 1** — version becomes `0.1.0` if we tag, but the "tag every phase" line in the roadmap is still aspirational; tagging cadence will be settled when (or if) we cut the first pre-release.
