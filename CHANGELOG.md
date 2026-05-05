# Changelog

All notable changes to gridcalc are documented in this file. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project adheres to [Semantic Versioning](https://semver.org/).

> **Roadmap renumber, 2026-05-04.** A new **Phase 14 — Finite volume
> method (FVM)** was inserted into `specs/roadmap.md` between (the
> just-merged) Phase 13 (Vector and tensor fields) and the
> previously-numbered Phase 14 (Functional evaluation for vector/tensor
> data, now **Phase 15**). All subsequent phases shifted up by one:
> `Phase N → Phase N+1` for `N ∈ [14, 22]`. The mechanical rename was
> applied to `specs/roadmap.md`, `specs/STATUS.md`, the post-Phase-13
> doc chapters' forward-references, and a single comment in
> `test/integrate_test.cpp`. Historical CHANGELOG entries below
> (0.11.0 / 0.12.0 / 0.13.0) are **not** renumbered — they are
> version-block-scoped historical records, and the references to
> "Phase 14" in those blocks meant the originally-numbered Phase 14
> (= today's Phase 15). Future blocks reference the post-renumber
> numbering directly.

## 0.14.2 — Fix: silent Grid-mismatch on RHS return in time integrators

### Fixed

- **Silent `Grid` mismatch on RHS return in `solve::integrate(...)`.**
  Both integrator overloads (`ExplicitEuler` and `RK4`) call the
  user-supplied `rhs(psi)` callable each stage and then loop
  `n = 0..psi.getSize() - 1` over the returned field's `data()`
  pointer. The compile-time check
  `std::is_invocable_r_v<Field<double>, Rhs&, const Field<double>&>`
  only constrains the *return type*, not the runtime `Grid` carried by
  the returned field. If the RHS allocates against a different `Grid`,
  three failure modes silently occurred — out-of-bounds read on a
  smaller field, layout drift across same-total-size / different-
  per-axis-extent fields, and `1/h`-magnitude drift across same-shape /
  different-cell-size fields. None tripped any existing test or warning.

  Both integrators (and the shared `solve::detail::axpyFresh` stage
  builder) now validate that `rhs(...).getGrid()` bit-matches
  `psi.getGrid()` immediately after each call and throw
  `std::invalid_argument` if not. Error messages name the offending
  RK4 stage (`k1` / `k2` / `k3` / `k4`) or the ExplicitEuler `rhs(psi)`
  call site, and print both the expected and actual
  `(Nx, Ny, Nz, hx, hy, hz)` tuples. The bug was dormant in this
  codebase — every RHS produced by the project's own `diff::*` /
  `fvm::*` operators is grid-preserving by construction — but the
  mission target ("production / industrial") includes user code that
  may compose RHS callables holding captured grids; surfacing the
  contract loudly catches the mistake at the integrator boundary
  rather than as a corrupted state vector mid-run.

### Added

- **`gridcalc::core::Grid::operator==` and `operator!=`** —
  bit-exact componentwise equality on `(Nx, Ny, Nz)` and the three
  components of `cell_size`. Used by the new precondition; documented
  as "bit-exact, by design" so a future reader does not relax it to a
  tolerance-based comparison (this codebase has no source of `Grid`
  numerical drift; any non-bit-exact mismatch is by definition a
  different `Grid` object). Carries `\since 0.14.2`.

### Tests

- New `test/grid_equality_test.cpp` (10 tests): identical-grid
  equality, reflexive / copy equality, per-axis extent mismatch
  (Nx, Ny, Nz), per-component cell-size mismatch (hx, hy, hz), and
  the bit-exact contract documented via `0.1 + 0.2 != 0.3`.
- New `test/integrator_grid_mismatch_test.cpp` (14 tests): three
  flavors of RHS-grid mismatch (different total cell count, same
  total / different per-axis extents, same shape / different cell
  size) on each of `RK4` and `ExplicitEuler`; direct coverage of
  `detail::axpyFresh` for matched and mismatched grids; per-stage
  label sweep on RK4 (`k1` / `k2` / `k3` / `k4`) verifying the helper
  names the offending stage in its message; and happy-path smoke
  tests confirming a `diff::laplacian`-based RHS integrates without
  throwing through both integrator overloads.
- Test totals: `clang-debug` is now 316 (292 prior + 24 new);
  `clang-debug-nofft` is now 240 (216 prior + 24 new).

### Internal

- New header `include/gridcalc/solve/detail/RhsGridCheck.hpp` with
  the free function
  `gridcalc::solve::detail::requireSameGrid(expected, actual, context)`.
  Lives under `solve/detail/` so it is excluded from the public-API
  surface and the Phase 10 `\since`-tag lint regex; called from
  `solve::detail::axpyFresh` and from each per-stage check inside
  `integrate(... RK4)` / `integrate(... ExplicitEuler)`. Carries
  Doxygen with `\since 0.14.2`.
- `axpyFresh`, `integrate(... RK4)`, `integrate(... ExplicitEuler)`,
  and the Phase 5 thin wrapper `explicitEuler` each gain a `\throws`
  line on the existing Doxygen block; `\since` tags note the
  precondition addition (e.g.,
  `\since 0.6.0 (function); 0.14.2 (RHS-grid precondition).`).

## 0.14.1 — Fix: silent stencil aliasing on small periodic axes

### Fixed

- **Silent stencil aliasing on small periodic axes.** Centered
  finite-difference stencils of radius `r` (Phase 7's
  `Coefficients<4>`, `FirstDerivative<4>`; Phase 8's
  `ThirdDerivative<{2, 4}>`, `FourthDerivative<{2, 4}>`) reach offsets
  `[-r, r]` along each active axis. On a periodic axis with extent
  `2 <= N <= 2 * r`, two or more offsets aliased to the same grid
  sample after the periodic wrap; the aliased weights summed at the
  duplicated index, silently degrading the operator away from its
  tabulated truncation order. The Phase 9 FD–FFT cross-check fixture's
  grid sweep `{16, 32, 64, 128}` was well above every threshold and
  did not catch this. Concrete failure case from the bug report:
  `diff::laplacian<4>(field)` on a grid with `Nx = 4` collapses the
  five-point `O(h^4)` stencil to a four-point operator with
  `{4/3, -8/3, 4/3, -1/6}` weights — published convergence order no
  longer holds.

  Every operator in `gridcalc::diff::*` and `gridcalc::fvm::*` that
  reads off-axis samples through the periodic index policy now
  validates each axis's extent at the call entry and throws
  `std::invalid_argument` if any axis falls in the partial-alias range
  `[2, 2 * radius]`. Affected operators: `diff::laplacian<{2, 4}>`,
  `diff::gradient<{2, 4}>` (scalar and Phase 13 vector overload),
  `diff::divergence<{2, 4}>` (scalar-output and Phase 13 tensor-output
  overload), `diff::biharmonic<{2, 4}>`, every Phase 8 d-prefix mixed
  partial at both orders (via the shared `mixedDerivative` site), and
  `fvm::cellLaplacian` (defense-in-depth — radius `1` at Phase 14, so
  currently a no-op for any `N >= 3`; the check is forward-compatible
  with any future higher-radius FVM variant).

  **Degenerate-axis carve-out.** `N = 1` is accepted regardless of
  stencil radius. At `N = 1` every offset wraps to index `0`, all
  reads return the single cell value, and the centered FD weight
  tables in this codebase sum to zero by construction (the n-th
  derivative of a constant is zero); the aliased result is therefore
  exactly `0` — the correct value for a degenerate single-cell
  periodic axis with no spatial variation. Practical consequence: the
  existing 1D-flavored `Grid g(4, 1, 1, ...)` test pattern used by
  Phase 14's FVM unit tests (and by user code with degenerate axes)
  continues to work without modification.

  Error messages name the offending axis (`x` / `y` / `z`), its `N`,
  the stencil's radius, and the required minimum.

### Tests

- New `test/stencil_axis_extent_test.cpp` (34 tests):
  - **Throw tests** at `N = 2 * radius` for each radius tier
    (`1, 2, 3`) across every affected operator.
  - **Per-axis label sweep** verifying the helper reports the right
    axis (`x`, `y`, `z`) in its error message.
  - **Smallest-valid-N boundary tests** at `N = 2 * radius + 1` for
    every Phase 7 / Phase 8 / Phase 13 / Phase 14 operator, asserting
    no-throw and finite output. The Phase 9 FD–FFT cross-check
    fixture's slope sweep is left at `{16, 32, 64, 128}` — extending
    it to `N = 5` is unsound (the truncation-error analysis breaks
    down at `h ≈ 1.257`); the no-throw + finiteness gate captures the
    actual silent-degradation failure mode without making slope
    claims at degenerate resolutions.
  - **Degenerate-axis carve-out tests** verifying that `N = 1` is
    accepted across `laplacian<4>`, `gradient<4>`, `biharmonic<4>`,
    `fvm::cellLaplacian`.
  - **Inactive-axis carve-out test** through the Phase 8 mixed-partial
    path (`mixedDerivative<3, 0, 0, 2>` on a grid with `Ny = Nz = 1`),
    confirming an axis with per-axis derivative count `0` accepts any
    positive `N`.
- Test totals: `clang-debug` is now 292 (258 prior + 34 new);
  `clang-debug-nofft` is now 216 (182 prior + 34 new).
- Phase 9 FD–FFT cross-check fixture is unchanged.

### Internal

- New header `include/gridcalc/diff/detail/PreconditionAxisExtent.hpp`
  with the free function
  `gridcalc::diff::detail::requireAxisExtent(axis_label, N, radius)`.
  Lives under `detail/` so it is excluded from the public-API surface
  and the Phase 10 `\since`-tag lint regex; reused by
  `fvm::cellLaplacian` via include. Carries Doxygen with `\since 0.14.1`.
- Each touched operator's existing Doxygen block gains a `\throws`
  line referencing the precondition; `\since` tags note the precondition
  addition (e.g., `\since 0.1.0 (function); 0.7.0 (Order parameter); 0.14.1 (precondition).`).

## 0.14.0 — Phase 14: Finite volume method (FVM)

### Added

- **`gridcalc::fvm` namespace and `include/gridcalc/fvm/CellLaplacian.hpp`**
  — new public header. Exposes
  `fvm::cellLaplacian(psi, D)` returning `∇·(D ∇ψ)` via cell-flux
  discretization with face-centered **harmonic-mean** D-averaging
  (`D_face = 2 D_L D_R / (D_L + D_R)`). Strictly positive `D` required
  per the harmonic-mean contract. The constant-`D = c` case reduces to
  `c · diff::laplacian(psi)` to round-off, verified by a dedicated
  agreement-gate test.
- **Heterogeneous-D `solve::diffuse` overload** —
  `solve::diffuse(psi, const Field<double>& D, dt, n_steps, tag)`.
  Reuses the existing `solve::ExplicitEuler{}` / `solve::RK4{}`
  integrator tags. Pre-loop pass over `D` validates the
  `D > 0` contract and computes `D_max` for the generalized CFL bound
  in a single sweep; throws `std::invalid_argument` with a flat-index
  pointer on contract violation. CFL bound generalizes to
  `D_max · dt · sum_a (1/h_a²) ≤ Tag::diffusionCFLLimit`. The
  constant-`D` Phase 5 / Phase 6 overload is unchanged.
- **Textbook FVM demo** — `examples/heterogeneous_diffusion.cpp`
  solving heterogeneous-`D` diffusion of a Gaussian peak through a
  smooth `D(x) = D0 + ε cos(2π x / L)` field. CLI options for grid
  size, time step (auto-CFL when `--dt 0`), step count, snapshot
  cadence, D parameters, IC sigma, integrator (`euler` | `rk4`), and
  output dir. VTK `STRUCTURED_POINTS` snapshots written via the Phase
  12 writer. Per-snapshot diagnostics on stdout: mass (constant to
  round-off), L²-norm-squared (strictly decreasing), psi_min (≥ 0),
  psi_max (monotonically non-increasing).
- **Demo helpers** — `examples/common/heterogeneous_diffusion.hpp`
  shared between the demo binary and its CI gate. Exposes
  `makeGaussianIc`, `makeHeterogeneousD` (with `|ε| < D0` validation),
  `computeMass`, `computeL2Squared`, `computeMin`, `computeMax`.
- **15 new tests** across four files:
  - `test/fvm_cell_laplacian_test.cpp` (6 tests): constant-D agreement
    with the FD path; harmonic-mean face flux on a 4-cell D-jump
    test; order-2 convergence on a manufactured `D(x) = 1 + 0.5 cos(x)`,
    `ψ(x) = sin(x)` pair; mass conservation on the periodic sum.
  - `test/diffusion_test.cpp` (4 tests): heterogeneous-D path agrees
    with constant-D path on a uniform `D` field; CFL violation throws;
    non-positive `D` throws; mass conserved over many integrator
    steps.
  - `test/fvm_diffusion_acceptance_test.cpp` (2 tests): order-2
    convergence on the analytical sin·sin·sin eigenfunction (uniform-D
    code path); qualitative properties on a genuinely varying `D(x)`
    (mass conservation, max|ψ| monotone-decrease, positivity
    preservation).
  - `test/heterogeneous_diffusion_test.cpp` (3 tests): demo CI gate
    with static-cached snapshot fixture (mirrors Phase 12's CH test
    pattern); ~1.6 s on Apple Clang debug.
- **User Guide chapter 14** — *Finite volume method* —
  `docs/user-guide/chapters/14-finite-volume-method.tex`. Walks
  motivation (FVM vs. FD), the cell-flux discretization, harmonic-mean
  averaging, the heterogeneous-D `solve::diffuse` API, and the demo
  walkthrough with three committed PNG snapshots showing anisotropic
  Gaussian spreading through the heterogeneous medium.
- **Developer Note chapter 13** — *Finite volume method* —
  `docs/developer-note/chapters/13-finite-volume-method.tex` with the
  standing five-section structure. Theory: conservation-law form of
  the PDE; Patankar-style derivation of the harmonic mean from the
  steady-state-flux match across a discontinuity. Math derivation:
  cell-flux discretization, order-2 truncation, generalized CFL bound.
  Algorithm: file layout, `cellLaplacian` loop, heterogeneous-D
  solver. Design decisions: the four AskUserQuestion answers + the
  manufactured-solution trap from Group 6. References: Patankar 1980,
  LeVeque 2007, Hundsdorfer & Verwer 2003, VTK file-format
  reference (all with permanent identifiers / DOIs / ISBNs).
- **Roadmap renumber.** New Phase 14 entry inserted into
  `specs/roadmap.md`; existing Phases 14–22 renumbered to 15–23.
  v1.0 release becomes Phase 23.
- **`scripts/render_diffusion_snapshots.py`** — matplotlib helper
  generating `z`-midplane PNGs from the demo's VTK snapshots
  (mirrors `render_ch_snapshots.py` from Phase 12 with a
  sequential-colormap default for non-negative `ψ`).

### Changed

- **Version bumped to `0.14.0`.**
- `docs/user-guide/main.tex` — adds
  `\input{chapters/14-finite-volume-method.tex}`.
- `docs/developer-note/main.tex` — adds
  `\input{chapters/13-finite-volume-method.tex}`.
- `examples/CMakeLists.txt` — new `heterogeneous_diffusion` target
  alongside `cahn_hilliard` (same gating, same warning set, picks up
  `examples/common/` via the shared `target_include_directories`).

### Notes

- **Hybrid incorporation, not replacement.** Phase 1's
  `diff::laplacian`, Phase 2's `gradient` / `divergence`, and Phase
  5 / Phase 6's constant-D `solve::diffuse` are unchanged.
- **Phase 9's FD–FFT cross-check fixture is unaffected** — Phase 14
  adds no new FD operators.
- **Periodic only.** Phase 19 (renumbered) handles non-periodic BCs
  and now leverages the FVM foundation.
- **No higher-order FVM, no IMEX.** 2nd-order central-difference face
  fluxes only; explicit time stepping only.

### Tests

- 258 tests pass on `clang-debug` (243 prior + 15 new spanning four
  files).
- 182 expected on `clang-debug-nofft` (167 prior + 15 new — none of
  the new tests touch the spectral path).

## 0.13.0 — Phase 13: Vector and tensor fields

### Added

- **`gridcalc::core::Mat3d`** alias for `Eigen::Matrix3d` in
  `core/EigenAliases.hpp`. Standard Eigen module — the unsupported
  tensor module remains deferred. Doxygen documents the index
  convention adopted by Phase 13's rank-2 gradient
  (`M(i, j) = ∂_j v_i`, so the trace gives the divergence) and the
  Eigen contract that the default constructor leaves coefficients
  uninitialized.
- **Rank-raising gradient overload** —
  `diff::gradient(Field<Vec3d>) -> Field<Mat3d>`. Returns the
  velocity-gradient tensor `M(i, j) = ∂_j v_i` at every grid point.
  Built as a sibling overload to the Phase 2 scalar gradient with the
  Phase 7 `Order` template parameter preserved (default `2`,
  `4` available).
- **Rank-lowering divergence overload** —
  `diff::divergence(Field<Mat3d>) -> Field<Vec3d>`. Returns
  `(div M)_i = ∂_j M(i, j)` (sum over the column / axis index, matching
  the rank-2 gradient convention so that
  `divergence(gradient(v))` is the vector Laplacian). Same
  `Order` template parameter.
- **`tensor/Contraction.hpp`** — three pointwise primitives in the new
  `gridcalc::tensor` namespace: `trace(Field<Mat3d>)`,
  `singleContract(Field<Mat3d>, Field<Mat3d>) -> Field<Mat3d>` (matrix
  product), and `doubleContract(Field<Mat3d>, Field<Mat3d>) -> Field<double>`
  (full contraction `A:B = A(i,j)B(i,j)`). All three are
  fully-materialized; expression templates revisited in Phase 14.
- **Vector-identity acceptance tests** (`test/vector_tensor_test.cpp`,
  16 tests). Covers per-operator unit checks (linear-field exactness,
  index convention, `Order=4` improvement) and three vector
  identities: `HessianIsSymmetric` (the rank-raising restatement of
  the roadmap's `∇×∇φ = 0` example), `TraceOfGradientEqualsDivergence`,
  and `DivergenceOfGradientEqualsLaplacian` (compared against the
  *analytical* vector Laplacian per the standing
  divergence-of-gradient stride-2-stencil note from Phase 2).
- **`Field<Mat3d>` instantiation coverage** (`field_test.cpp`, four
  new tests): broadcast constructor, zero broadcast, write-then-read
  round-trip, i-fastest layout preserved at the `Mat3d` granularity.
- **User Guide chapter 13** — *Vector and tensor fields* —
  `docs/user-guide/chapters/13-vector-tensor-fields.tex`. Walks the
  new types and overloads, the index convention, the contraction
  primitives, and the three identity tests. The legacy
  `13-diamond-lattice.tex` placeholder (a Phase 11 renumber-residue
  for Phase 15's diamond-lattice chapter) moves to
  `15-diamond-lattice.tex`.
- **Developer Note chapter 12** — *Vector and tensor fields* —
  `docs/developer-note/chapters/12-vector-tensor-fields.tex` with the
  standing five-section structure (Theory · Math derivation ·
  Algorithm · Design decisions · References). Theory: continuum
  rank-raising gradients and Schwarz's theorem. Math derivation:
  proof of discrete partial-difference commutativity on a periodic
  uniform grid, the trace=divergence pointwise-to-round-off
  argument, and the stride-2 second-derivative analysis for
  `divergence(gradient(v))`. Algorithm: file layout, loop structure
  for the rank-2 gradient, and the contraction-primitive
  implementation. Design decisions: the four AskUserQuestion
  answers. References: Gurtin 1981, Holzapfel 2000, LeVeque 2007,
  Eigen project documentation (all with permanent identifiers /
  ISBNs).

### Changed

- **Version bumped to `0.13.0`.**
- `docs/user-guide/main.tex` and `docs/developer-note/main.tex` —
  add the new chapter `\input{...}` lines; legacy
  `13-diamond-lattice.tex` reference replaced with
  `15-diamond-lattice.tex` to free chapter slot 13 for the Phase 13
  content. (Phase 14 will land `14-tensor-functional.tex` between
  them.)

### Notes

- **No new public namespace pollution.** All new symbols live in
  existing namespaces (`gridcalc::core`, `gridcalc::diff`) plus the
  one new `gridcalc::tensor` namespace introduced by the new header.
  The `\since` lint passes — every new public symbol carries a
  `\since 0.13.0` tag.
- **No new FD operators in the Phase 9 cross-check fixture.** The
  Phase 9 fixture parameterizes over scalar-input FD operators; the
  rank-raising variants are not (yet) in scope. May be addressed in a
  Phase 14-or-later body of work.
- **Phase 12's CH demo unchanged.** The new vector/tensor surface is
  not consumed by Phase 12; Phase 14's elastic-energy worked example
  is the first downstream user.

### Tests

- 243 tests pass on `clang-debug` (223 prior + 20 new — 4
  `FieldMat3dTest` + 13 in `vector_tensor_test.cpp` + 3 identity
  tests in the same file).
- 167 expected on `clang-debug-nofft` (147 prior + 20 new — none of
  the new tests use the spectral path).

## 0.12.0 — Phase 12: Cahn–Hilliard example end-to-end

### Added

- **`examples/` subtree, gated by a new `GRIDCALC_BUILD_EXAMPLES`
  CMake option** (default `OFF`, mirroring Phase 10's
  `GRIDCALC_BUILD_DOCS`). The first example, `cahn_hilliard`, solves
  the conservative phase-field equation
  $\partial_t\psi = M\nabla^2(-\psi + \psi^3 - \kappa\nabla^2\psi)$
  on a periodic 3D grid using explicit RK4 (Phase 6
  `solve::integrate(..., RK4{})`) and the 2nd-order finite-difference
  Laplacian (Phase 1 `diff::laplacian`). CLI options expose grid size,
  time step, integration length, snapshot cadence, $\kappa$ / $M$,
  RNG seed, and output directory; the binary writes one VTK snapshot
  per `--snapshot-every` steps and prints free energy / gradient
  energy / mean / `max|psi|` to stdout. CI builds the example under
  the project's strict warning set on every PR (`-DGRIDCALC_BUILD_EXAMPLES=ON`)
  but does not run it; the runtime gate is the new acceptance test.
- **`examples/common/cahn_hilliard.hpp`** — header-only RHS and energy
  diagnostics shared between the demo binary and the acceptance test.
  Exposes `computeRhs`, `computeFreeEnergy`, `computeGradientEnergy`,
  `computeMean`, and `makeRandomIc` under
  `gridcalc::examples::cahn_hilliard`. The energy helpers route
  through the Phase 4 / Phase 11 `func::evaluate` API to exercise the
  functional surface on a real workload.
- **`examples/common/vtk_writer.hpp`** — hand-written legacy ASCII
  VTK 3.0 `STRUCTURED_POINTS` writer for `Field<double>`. ~30 lines,
  no new build dependency, opens directly in ParaView and VisIt.
  `gridcalc::examples::writeVtkStructuredPoints(field, path, name)`
  emits the `DIMENSIONS` / `ORIGIN` / `SPACING` block followed by a
  single `SCALARS <name> double 1` array; the i-fastest storage layout
  matches VTK's expected order for `STRUCTURED_POINTS` so the payload
  is a flat read of `field.data()`.
- **Acceptance test** (`test/cahn_hilliard_test.cpp`, six tests).
  Three quick correctness checks (`ZeroRhsForConstantField`,
  `RhsIsConservative`, `SeedReproducibility`) plus three properties of
  a single ~30 s simulation shared across the three TEST_F bodies via
  a static-cached lambda: `LyapunovDecay` (F(t) is monotonically
  non-increasing), `MonotoneCoarsening` (gradient energy
  $\int|\nabla\psi|^2\,dV$ drops at least 3% from t=30 (peak) to t=50
  (post-coarsening); locally observed 4.83%–7.85% across seeds {1, 2,
  7, 42, 1234}), `MassConservation`
  (`|<psi>(t_final) - <psi>(0)| < 1e-10`). Discharges the roadmap
  acceptance bar "coarsening dynamics qualitatively match published
  CH benchmarks."
- **VTK writer round-trip test** (`test/vtk_writer_test.cpp`, three
  tests): header-shape, payload exact-match in i-fastest order, and
  scalars-line field-name round-trip.
- **`scripts/render_ch_snapshots.py`** — matplotlib helper that parses
  the legacy ASCII VTK snapshots and renders 2D z-midplane PNGs for
  the User Guide chapter. Headless ("Agg" backend); used to generate
  the three committed coarsening figures from a 64³ release-build run
  at the demo's default parameters.
- **User Guide chapter 12** — *Cahn–Hilliard demo* —
  `docs/user-guide/chapters/12-cahn-hilliard.tex` (replaces the Phase
  11 placeholder). Walks the equation → code → results sequence with
  the closed-form variational derivative, the $\Delta x^4$ explicit
  RK4 stability budget, a code-walkthrough citing `computeRhs` and the
  `func::evaluate`-routed energy diagnostics, run instructions, and a
  three-figure coarsening sequence (early / mid / late at t=20, 40,
  80).
- **Developer Note chapter 11** — *Cahn–Hilliard demo* —
  `docs/developer-note/chapters/11-cahn-hilliard.tex` with the
  standing five-section structure (Theory · Math derivation ·
  Algorithm · Design decisions · References). Theory: Cahn–Hilliard
  1958 I+II, the Lyapunov property, Bray's classical-coarsening
  scaling. Math derivation: Euler–Lagrange δF/δψ, linearized
  dispersion σ(k) = M k²(1 − κk²), discrete worst-case eigenvalue
  feeding the RK4 Δt budget. Algorithm: file layout, build wiring,
  per-step allocation budget (~20 fields / step at the default
  resolution), VTK writer format, acceptance-test instrumentation.
  Design decisions: the four AskUserQuestion answers + helpers-under-
  examples decision. References: Cahn–Hilliard 1958 I+II, Bray 1994,
  Provatas–Elder 2010, LeVeque 2007, Hairer–Nørsett–Wanner 1993, and
  the VTK file-format reference (all with permanent
  identifiers / DOIs).

### Changed

- **Version bumped to `0.12.0`.**
- **CI configures with `-DGRIDCALC_BUILD_EXAMPLES=ON`** on both Ubuntu
  jobs so the example is type-checked under the project's strict
  warning set on every PR.
- **`test/CMakeLists.txt`** adds `${CMAKE_SOURCE_DIR}/examples` to
  `gridcalc_tests`'s include path so the helper headers under
  `examples/common/` are reachable from `test/` as
  `#include <common/...>`.

### Notes

- **No new public symbols in `include/gridcalc/`.** Phase 12 consumes
  Phases 1–8 and adds files only under `examples/`, `test/`, `docs/`,
  and `scripts/`. The Phase 10 `\since` lint regex covers
  `include/gridcalc/` only and is unaffected. Phase 9's FD–FFT
  cross-check fixture is also unaffected (no new FD operators).
- **No FFT-diagonalized / IMEX time stepping.** The "FFT verification
  only" rule (mission.md) and the deferral of implicit schemes to
  Phase 19 keep that work out of this phase. The
  $\Delta t \sim h^4 / (M\kappa)$ explicit-RK4 stability cost is
  documented in both doc chapters.

### Tests

- 223 tests pass on `clang-debug` (214 prior + 6 new in
  `cahn_hilliard_test.cpp` + 3 new in `vtk_writer_test.cpp`).
- 147 tests pass on `clang-debug-nofft` (138 prior + 9 new; none of the
  new tests use the spectral path so the FFT-off build picks up all
  nine).
- The CH gate runs in ~30 s on Apple Clang debug; well under the
  60 s budget.

## 0.11.0 — Phase 11: Functional with up-to-4th-order gradients

### Added

- **`gridcalc::func::evaluate` 4-argument arity.** New SFINAE-detected
  branch at the head of the dispatch ladder accepts callables of the
  form `f(double psi, const Vec3d& grad_psi, double lap_psi, double
  biharm_psi)`, where `biharm_psi` is the contracted scalar
  $\nabla^4\psi$ as returned by `diff::biharmonic` (Phase 8). When this
  arity matches, `evaluate` materializes `gradient(psi)`,
  `laplacian(psi)`, and `biharmonic(psi)` in addition to the integrand
  field, then forms the integrand pointwise and reduces with
  `func::integrate`. Total supported arities: 1, 2, 3, 4.
- **PFC functional test suite** (`test/pfc_functional_test.cpp`,
  five tests): `HandComputedAtN64` (roadmap acceptance — recovers
  $F_{\text{ref}} = 2.05546875\,\pi^3$ within 2% relative at `N=64`
  for `psi = sin(x)sin(y)sin(z)`, `q0=1`, `eps=0.1`),
  `SecondOrderConvergence` (`O(h²)` slope on `N ∈ {16, 32, 64}`),
  `PsiOnlyArityStillWorks` and `GinzburgLandauStillWorks`
  (regression checks confirming the 1-arg and 3-arg dispatch branches
  still resolve correctly under the longest-arity-wins ordering), and
  `SpectralCrossCheck` (computes the same PFC integrand using
  Phase 9's `spectral::laplacian` / `spectral::biharmonic` and asserts
  agreement with the analytical reference to better than `1e-12`
  relative; `#ifdef`-guarded by `GRIDCALC_USE_FFT`).
- **User Guide chapter 11** — *Higher-order functionals* —
  `docs/user-guide/chapters/11-higher-order-functional.tex`. Walks
  the new arity, the PFC worked example end-to-end (parameters →
  integrand → analytical reference → calling code), and the spectral
  cross-check pattern. Renumbers the placeholder chapters: previous
  chapter 11 (Cahn–Hilliard) becomes chapter 12; previous chapter 12
  (diamond-lattice) becomes chapter 13.
- **Developer Note chapter 10** — *Higher-order functionals* —
  `docs/developer-note/chapters/10-higher-order-functional.tex`. The
  five-section structure: Theory (PFC functional and the
  Swift–Hohenberg lineage of $(q_0^2+\nabla^2)^2$); Math derivation
  (closed-form $F_{\text{ref}}$, leading-order discrete error
  $\sim h^2$); Algorithm (dispatch ladder extension, eager
  materialization, spectral cross-check assembly); Design decisions
  (contracted-scalar-only argument shape, deferral of contraction
  expression templates to Phase 14, hand-computed + spectral
  verification); References (Elder & Grant 2004, Provatas & Elder
  2010, Boyd 2001, Swift & Hohenberg 1977).

### Changed

- `func::evaluate`'s `\since` tag now reads `0.4.0 (function); 0.11.0
  (4-arg arity)`. The `\file`-level brief is updated to reflect the
  $\nabla^4\psi$ argument; the `static_assert` error message
  enumerates all four supported arities.

## 0.10.0 — Phase 10: Documentation infrastructure

### Added

- **Doxygen + LaTeX documentation toolchain.** New `docs/Doxyfile.in`
  configured for the public headers (excluding `detail/`), with
  Graphviz `dot` enabled for class / include graphs. New `docs/CMakeLists.txt`
  declares three custom targets — `gridcalc-docs-doxygen`,
  `gridcalc-docs-userguide`, `gridcalc-docs-developernote` — and an
  umbrella `gridcalc-docs-all`. New build flag `GRIDCALC_BUILD_DOCS`
  (default `OFF` so contributors without TeX Live aren't blocked).
- **User Guide LaTeX skeleton** (`memoir` class) under
  `docs/user-guide/`. New "Hello grid" tutorial chapter walks through
  the operators delivered through Phase 9 with one continuous example.
  Per-phase chapters 02–10 carry the absorbed-from-markdown reference
  content; chapters 11–12 are placeholder stubs reserved for Phases 12
  (Cahn–Hilliard) and 16 (diamond-lattice).
- **Developer Note LaTeX skeleton** (`book` class) under
  `docs/developer-note/`. Chapters 01–09 carry the absorbed
  developer-note content (five-section structure preserved). The back
  matter `\input{}`s the Doxygen-generated content extracted from
  `refman.tex` via `docs/extract-api-reference.sh` so the auto-built
  API reference for every public Phase 0–9 symbol lives in the same
  PDF.
- **Per-phase markdown notes absorbed into LaTeX chapters.** The 18
  files under `docs/{user-guide,developer-note}/notes/phase-N-*.md`
  are pandoc-converted into LaTeX chapters; the `notes/` directories
  are emptied. From Phase 11 onward, doc deliverables go straight to
  the LaTeX sources.
- **`scripts/build-docs.sh`** — driver that runs the full docs pipeline
  via the CMake targets above.
- **`scripts/check-since.py`** — `\since`-tag lint over public headers.
  Walks `include/gridcalc/**/*.hpp` (excluding `detail/`), identifies
  public-API declarations, and verifies each preceding Doxygen block
  contains an explicit `\since` line. Exits non-zero on any miss.
- **CI `docs-build` job** installs TeX Live + Doxygen + Graphviz and
  produces three artifacts in the `gridcalc-docs-pdfs` bundle (Doxygen
  output tarball, User Guide PDF, Developer Note PDF) on every PR.
- **CI `docs-lint` job** runs Doxygen with `WARN_IF_UNDOCUMENTED=YES`
  and parses the warnings log (allowing benign math-mode warnings
  only); also runs `scripts/check-since.py`. Becomes a permanent CI
  gate from Phase 10 forward.

### Changed

- `specs/tech-stack.md` Documentation section updated: pandoc removed
  from the required toolchain (was only used by the lightweight
  pre-Phase-10 render); the `\since` policy is now CI-enforced.
- `scripts/get-docs.sh` `--local` mode now drives the new full-LaTeX
  pipeline via `scripts/build-docs.sh` instead of the retired pandoc
  render. The `gridcalc-docs-pdfs` artifact name is unchanged.
- `include/gridcalc/solve/{ExplicitEuler,RK4}.hpp` Doxygen comments
  switched from `$\dot\psi = ...$` to `$\partial_t\psi = ...$`. The
  earlier form tripped Doxygen's `\dot` graph-block parser; the math
  is equivalent.

### Removed

- **`scripts/render-docs.sh`** — retired alongside the CI
  `Render Doc PDFs` job in favor of the full LaTeX skeleton.
- **`docs/{user-guide,developer-note}/notes/phase-N-*.md`** (18 files)
  — absorbed into LaTeX chapters.

### Documentation

- `docs/README.md` documents the local build commands for all three
  artifacts, the required local toolchain, and the `scripts/get-docs.sh --ci`
  fallback.

## 0.9.0 — Phase 9: Spectral verification harness

### Added

- **Vendored PocketFFT** under `extern/pocketfft/pocketfft_hdronly.h`
  (BSD-3-Clause). Revision recorded in `extern/pocketfft/VERSION.txt`.
  Vendored rather than `FetchContent`-fetched (decided in the Phase 9
  spec round): build determinism without network access matters for the
  project's "production / industrial" target.
- New CMake option `GRIDCALC_USE_FFT` (default `ON`). When `ON`, the
  `gridcalc::spectral` namespace is built and the FD–FFT cross-check
  fixture is included in the test suite. When `OFF`, all spectral
  headers `#ifdef`-out and the FD code path is unaffected.
- New namespace `gridcalc::spectral` (header-only):
  - `spectral::Wavenumbers` — `kxRfft`, `kyFull`, `kzFull` builders for
    per-axis $k$-space arrays following the `numpy.fft`-style signed
    convention (rfft half-spectrum on the x-axis).
  - `spectral::Fft` — `forwardR2C` and `backwardC2R` thin wrappers over
    PocketFFT's r2c / c2r 3D transforms. Column-major byte strides
    matching `Field<double>`'s i-fastest layout; `axes = {2, 1, 0}` so
    the rfft compresses axis 0 (the x-axis).
  - `spectral::Derivatives`:
    - `spectral::partial<int Nx, int Ny, int Nz>(field) -> Field<double>`
      — generic multi-index spectral partial. Multiplies the spectrum
      by $(ik_x)^{N_x}(ik_y)^{N_y}(ik_z)^{N_z}$ and inverse-transforms;
      Nyquist mode zeroed on any axis with an odd derivative count.
    - `spectral::laplacian(field)` — fused $-|\mathbf{k}|^2$ multiplication.
    - `spectral::biharmonic(field)` — fused $|\mathbf{k}|^4$ multiplication.
    - `spectral::gradient(field) -> Field<Vec3d>` — one forward + three
      backward FFTs.

### Tests

- **`SpectralBasicTest`** (5 tests) — FFT round-trip; named-wrapper
  smoke tests for `laplacian`, `biharmonic`, `gradient` on
  $\sin(\mathbf{k}\cdot\mathbf{x})$ matching their analytical
  $-|\mathbf{k}|^2\psi$, $|\mathbf{k}|^4\psi$, $\mathbf{k}\cos(\mathbf{k}\cdot\mathbf{x})$
  references; Nyquist-zeroing audit on a one-mode-at-Nyquist input.
- **`FdFftCrossCheck`** (70 tests) — the permanent CI gate. For every
  Phase 1–8 FD operator at every advertised accuracy order, asserts the
  log-log slope of the FD–FFT discrepancy on $\psi = \sin(x + y + z)$
  over $N \in \{16, 32, 64, 128\}$ lies in
  $[\text{Order} - 0.5, \text{Order} + 0.5]$. Coverage:
  Phase 1 `laplacian` × 2 orders, Phase 2 `gradient` / `divergence` × 2
  orders, Phase 8's 31 d-prefix partials × 2 orders + `biharmonic` × 2
  orders.
- All Phase 1–8 tests continue to pass without modification (verified
  with both `GRIDCALC_USE_FFT=ON` and `=OFF`).

### Documentation

- New User Guide note
  `docs/user-guide/notes/phase-9-spectral-verification.md` covering the
  hybrid API, verification-only positioning, build flag, and a worked
  CI-style FD–FFT assertion example.
- New Developer Note
  `docs/developer-note/notes/phase-9-spectral-verification.md` with the
  five-section structure (Theory · Math derivation · Algorithm · Design
  decisions · References), wavenumber + Nyquist conventions, three
  external references (Trefethen DOI, Frigo–Johnson DOI, Boyd permanent
  URL) plus the PocketFFT upstream URL.
- `specs/tech-stack.md` updated to record PocketFFT as **vendored**
  rather than `FetchContent`-fetched, with the BSD-3-Clause license note
  corrected (was "MIT").

## 0.8.0 — Phase 8: Higher-order spatial derivatives (3rd, 4th)

### Added

- `gridcalc::stencil::ThirdDerivative<Order>` — central-difference weight
  table for $\partial^3/\partial x^3$. Order 2: $\{-1/2, 1, 0, -1, 1/2\}$
  (radius 2, anti-symmetric). Order 4: $\{1/8, -1, 13/8, 0, -13/8, 1, -1/8\}$
  (radius 3, anti-symmetric).
- `gridcalc::stencil::FourthDerivative<Order>` — central-difference weight
  table for $\partial^4/\partial x^4$. Order 2: $\{1, -4, 6, -4, 1\}$
  (radius 2, symmetric). Order 4:
  $\{-1/6, 2, -13/2, 28/3, -13/2, 2, -1/6\}$ (radius 3, symmetric).
- `gridcalc::diff::detail::mixedDerivative<Nx, Ny, Nz, Order>` — internal
  tensor-product helper that computes the multi-axis partial
  $\partial^{Nx+Ny+Nz}\psi / \partial x^{Nx}\partial y^{Ny}\partial z^{Nz}$
  in a single fused pass via the outer product of three 1D weight tables.
- **31 named partial-derivative functions** across three new headers,
  each templated on `<int Order = 2>` (orders 2 and 4 specialized):
  - `gridcalc/diff/MixedPartial.hpp` — rank-2 partials: `dxx`, `dyy`,
    `dzz`, `dxy`, `dxz`, `dyz` (6 functions). Roadmap deliverable
    `diff/MixedPartial.hpp` shipped under the d-prefix naming
    convention.
  - `gridcalc/diff/ThirdOrder.hpp` — rank-3 partials (10 functions):
    `d3dx3`, `d3dy3`, `d3dz3`, `d3dx2dy`, `d3dxdy2`, `d3dx2dz`,
    `d3dxdz2`, `d3dy2dz`, `d3dydz2`, `d3dxdydz`.
  - `gridcalc/diff/FourthOrder.hpp` — rank-4 partials (15 functions):
    `d4dx4`, `d4dy4`, `d4dz4`, `d4dx3dy`, `d4dxdy3`, `d4dx3dz`,
    `d4dxdz3`, `d4dy3dz`, `d4dydz3`, `d4dx2dy2`, `d4dx2dz2`,
    `d4dy2dz2`, `d4dx2dydz`, `d4dxdy2dz`, `d4dxdydz2`.
- `gridcalc::diff::biharmonic<Order>(psi)` — $\nabla^4\psi$ on a periodic
  scalar field, implemented as a fused six-term direct stencil
  ($\partial_x^4 + \partial_y^4 + \partial_z^4 +
  2\sum_{a<b}\partial_a^2\partial_b^2$) — not via
  `laplacian(laplacian(psi))`. Roadmap deliverable
  `diff/Biharmonic.hpp`.
- `gridcalc::diff::d4<Order>(psi)` — alias forwarding to `biharmonic`,
  completing the d-prefix family at the contracted-rank level.

### Tests

- **Weight-table audits**: `EXPECT_DOUBLE_EQ` on every weight in
  `ThirdDerivative<{2, 4}>` and `FourthDerivative<{2, 4}>` (catches typos
  in the hard-coded values).
- **62 convergence sweeps** across all 31 named partials × `Order ∈ {2, 4}`
  on the trigonometric manufactured solution $\psi = \sin(x)\sin(y)\sin(z)$
  over the periodic $[0, 2\pi)^3$ box: log-log slope on $N \in \{16, 32, 64\}$
  in $[\text{Order} - 0.5, \text{Order} + 0.5]$. Satisfies the roadmap's
  acceptance.
- **Biharmonic convergence sweeps** (orders 2 and 4) and **four explicit
  $\nabla^4\psi = |\mathbf{k}|^4\psi$ acceptance tests** at $N = 64$ on
  $\mathbf{k} = (1, 1, 1)$ and $(1, 2, 1)$ — satisfies the roadmap's
  "biharmonic of $\sin(\mathbf{k}\cdot\mathbf{x})$ recovers $|\mathbf{k}|^4$"
  bar.
- **Alias-parity tests** verifying `d4<Order>(psi) == biharmonic<Order>(psi)`
  element-wise.
- All Phase 1–7 tests continue to pass without modification.

### Documentation

- New User Guide note
  `docs/user-guide/notes/phase-8-higher-order-derivatives.md` covering the
  full naming convention, the 31 named functions plus `biharmonic` / `d4`,
  the calling syntax with `Order`, and a worked $\nabla^4$ example with the
  $|\mathbf{k}|^4$ spectral expectation cross-check.
- New Developer Note
  `docs/developer-note/notes/phase-8-higher-order-derivatives.md` with the
  five-section structure (Theory · Math derivation · Algorithm · Design
  decisions · References), explicit Taylor-matching derivation of the four
  new 1D weight tables, the biharmonic six-term decomposition, and three
  external references (Fornberg 1988 DOI, LeVeque 2007 DOI, NIST DLMF
  permanent URL).
- `specs/roadmap.md` Phase 8 entry rewritten to match the expanded scope
  decided in the spec round.

## 0.7.0 — Phase 7: 4th-order accuracy stencils

### Added

- `gridcalc::stencil::Coefficients<4>` — 4th-order central-difference
  weights for $\partial^2/\partial x^2$:
  $\{-1/12, 4/3, -5/2, 4/3, -1/12\}$ (radius 2). Truncation error
  $O(h^4)$ vs `Coefficients<2>`'s $O(h^2)$.
- `gridcalc::stencil::FirstDerivative<4>` — 4th-order central-difference
  weights for $\partial/\partial x$:
  $\{1/12, -2/3, 0, 2/3, -1/12\}$ (radius 2). Truncation error
  $O(h^4)$ vs `FirstDerivative<2>`'s $O(h^2)$.

### Changed

- `gridcalc::diff::laplacian`, `gradient`, and `divergence` are now
  templated on an `Order` parameter, defaulting to `2`. All Phase 1–6
  callers (which write `laplacian(field)`, etc.) keep compiling
  unchanged via the C++17 default-template-argument rule. Phase 7
  callers can write `laplacian<4>(field)`, `gradient<4>(field)`, or
  `divergence<4>(vfield)` to opt into the 4th-order stencil.

### Tests

- Order-4 convergence sweep on each operator over `N ∈ {16, 32, 64}`:
  log-log slope in `[3.5, 4.5]`. Satisfies the roadmap's "log-log
  convergence slopes match advertised orders" acceptance.
- Order-2-vs-order-4 absolute-accuracy comparison on each operator at
  `N = 32`: order 4 is at least 10× more accurate than order 2.
- Weight-table verification for `Coefficients<4>` and
  `FirstDerivative<4>` (catches typos in the hard-coded values).

### Documentation

- New User Guide note `docs/user-guide/notes/phase-7-higher-order-stencils.md`
  with the API switch, cost / accuracy table, and a worked example.
- New Developer Note `docs/developer-note/notes/phase-7-higher-order-stencils.md`
  with the explicit Taylor / moment-matching derivation of both 4th-order
  weight tables, the truncation-error analysis, and four external
  references (Fornberg 1988, LeVeque, NIST DLMF, Trefethen).

## 0.6.0 — Phase 6: RK4 + generic time integrator interface

### Added

- `gridcalc::solve::ExplicitEuler` and `gridcalc::solve::RK4` empty
  tag structs, each carrying a `static constexpr double diffusionCFLLimit`
  (Euler `0.5`, RK4 `0.6963`).
- `gridcalc::solve::integrate(psi, rhs, dt, n_steps, IntegratorTag{})` —
  generic tag-dispatched time integrator. Two overloads, one per tag.
- Classic 4-stage Runge–Kutta implementation in
  `solve/RK4.hpp`; per-step allocation is `~6` fields. Phase 20 will
  optimize via persistent scratch.
- `solve/Integrator.hpp` — one-stop aggregator that pulls in both
  `ExplicitEuler.hpp` and `RK4.hpp`.

### Changed

- `gridcalc::solve::diffuse` is now templated on the integrator tag:
  `diffuse(psi, D, dt, n_steps, Integrator{})`. The default tag is
  `ExplicitEuler{}`, so Phase 5 callers (no integrator argument) keep
  working unchanged. The CFL check now reads
  `Integrator::diffusionCFLLimit`, so RK4 callers get the larger
  stability bound automatically.
- `gridcalc::solve::explicitEuler` is now a thin wrapper that forwards
  to `solve::integrate(..., ExplicitEuler{})`. Phase 5 callers and
  tests are unaffected.
- The Phase-5-specific `solve::detail::checkExplicitDiffusionCFL` is
  superseded by the templated `solve::detail::checkDiffusionCFL<Tag>`.

### Tests

- RK4 matches the analytical heat-equation decay within `1e-2` relative
  at `N=32, T=5.0` (10× tighter than Phase 5's Euler at the same
  parameters).
- Combined refinement (`N ∈ {16, 32, 64}` with each integrator's own
  CFL-limited dt): both Euler and RK4 show log-log slope in `[1.8, 2.5]`
  (combined order is spatial-dominated at CFL-limited dt).
- Pure linear-decay ODE comparison: RK4's error is at least 10× smaller
  than Euler's at the same dt (isolates time-error from spatial-error).
- **Pure linear-decay ODE convergence sweep over `dt`**: RK4 log-log
  slope in `[3.5, 5.0]`, ratio `>= 8` per halving. Satisfies the
  roadmap's `error scales as Δt^4` acceptance.
- RK4 succeeds in the gap between Euler's CFL (`0.5`) and RK4's
  (`0.6963`) where Euler throws.
- Bit-identical output between Phase 5's `solve::explicitEuler(...)`
  wrapper and `solve::integrate(..., ExplicitEuler{})`.

### Documentation

- New User Guide note `docs/user-guide/notes/phase-6-rk4-integrator.md`.
- New Developer Note `docs/developer-note/notes/phase-6-rk4-integrator.md`
  with the four-stage Taylor-match derivation, the RK4 stability-region
  intersection at `~-2.7853`, and four external references
  (Hairer–Nørsett–Wanner, Butcher, Wikipedia RK methods, LeVeque).

## 0.5.0 — Phase 5: explicit Euler diffusion solver

### Added

- `gridcalc::solve::explicitEuler<Rhs>(Field<double>&, Rhs&&, double dt, int n_steps)`
  — generic forward-Euler integrator; advances
  $\psi^{n+1} = \psi^n + dt \cdot \mathrm{rhs}(\psi^n)$ for an arbitrary
  RHS callable. Mutates `psi` in place. SFINAE-validates the callable
  signature; throws `std::invalid_argument` on negative `dt` or
  `n_steps`.
- `gridcalc::solve::diffuse(Field<double>&, double D, double dt, int n_steps)`
  — convenience driver for $\partial_t\psi = D\nabla^2\psi$. Internally
  calls `explicitEuler` with `diff::laplacian` as the RHS and `D * dt`
  as the effective step. Performs a von Neumann CFL check before
  integration:
  $D \cdot dt \cdot \sum_a (1/h_a^2) \le 1/2$.
  Throws `std::invalid_argument` on violation, with a message naming
  the offending values.
- New public namespace `gridcalc::solve` and directory
  `include/gridcalc/solve/`.

### Tests

- Trig eigenfunction `ψ_0 = sin(x)sin(y)sin(z)` decays to the analytical
  `exp(-3 D T) ψ_0` after 100 explicit-Euler steps; relative max-norm
  error `< 5e-2` at `N = 32`.
- 2nd-order convergence sweep on the same trig IC over `N ∈ {16, 32, 64}`
  with `dt ∝ h²`; log-log slope in `[1.8, 2.5]`.
- Gaussian sanity test: finite values, mass conserved (via
  `func::integrate`), peak max-norm decreased after integration.
- CFL violation throws.
- Generic `explicitEuler` with zero RHS leaves `psi` bit-identical.
- Negative `n_steps` throws on both entry points.

### Documentation

- New User Guide note `docs/user-guide/notes/phase-5-explicit-euler.md`
  with worked examples for `diffuse` and a custom-RHS use of
  `explicitEuler`.
- New Developer Note `docs/developer-note/notes/phase-5-explicit-euler.md`
  with the von Neumann derivation of the CFL bound and four external
  references (Strikwerda DOI, LeVeque DOI, Wikipedia heat-equation
  permanent URL, Numerical Recipes ISBN).

### Changed

- `specs/CLAUDE.md` step 4a: new "Commit-message rules" subsection
  forbidding `Co-Authored-By:` trailers in this repository.

## 0.4.0 — Phase 4: functional evaluation, callable API (scalar, periodic)

### Added

- `gridcalc::func::evaluate(const Field<double>&, F&&, Tag = {}) -> double`
  — discrete functional $F[\psi] = \int f(\psi, \nabla\psi, \nabla^2\psi)\, dV$.
  The callable `f` is auto-detected at compile time as one of three
  arities:
  - `f(double psi)` — only `ψ`; the gradient and Laplacian are skipped.
  - `f(double psi, const Vec3d& grad)` — `ψ` and `∇ψ`; the Laplacian is
    skipped.
  - `f(double psi, const Vec3d& grad, double lap)` — full triplet.

  The reduction tag (`Pairwise` or `Kahan`) is forwarded to
  `func::integrate`. Implementation is eager: required derivatives are
  materialized via `diff::gradient` / `diff::laplacian` and a
  per-grid-point integrand `Field<double>` is reduced.

### Tests

- Ginzburg–Landau hand-computed reference at `N = 64` matches
  $F = (507/64)\pi^3 \approx 245.62$ within `1e-2` relative.
- 2nd-order convergence sweep on the GL functional over `N ∈ {16, 32, 64}`.
- Per-arity dispatch tests: `f(ψ) = ψ²`, `f(ψ, ∇ψ) = (1/2)|∇ψ|²`,
  `f(ψ, _, ∇²ψ) = ψ ∇²ψ`, all matching their analytical values.
- Pairwise and Kahan reducers agree on the GL setup within `1e-13` relative.

### Documentation

- New User Guide note `docs/user-guide/notes/phase-4-functional-evaluation.md`
  walking through the Ginzburg–Landau worked example end-to-end (the
  "tutorial doc page" deliverable from `roadmap.md` Phase 4).
- New Developer Note `docs/developer-note/notes/phase-4-functional-evaluation.md`
  with the five-section structure: theory of functional convergence,
  full closed-form derivation of the GL reference value, arity-dispatch
  algorithm, design decisions, and five external references.

## 0.3.0 — Phase 3: domain integration (∫ over the grid)

### Added

- `gridcalc::func::integrate(const Field<double>&, Tag = {}) -> double` —
  discrete periodic-midpoint integral
  $I_h = h_x h_y h_z \sum_{i,j,k} f_{ijk}$. Two reduction tags:
  - `gridcalc::func::Pairwise` (default): pairwise recursive summation,
    rounding error $O(\varepsilon \log n)$, 64-element base case.
  - `gridcalc::func::Kahan`: Neumaier-style compensated summation,
    rounding error independent of $n$ to first order.
- New public namespace `gridcalc::func` and corresponding directory
  `include/gridcalc/func/`.

### Tests

- Unit-field-recovers-domain-volume on an anisotropic non-cubic Grid.
- $\int \sin(kx)\, dV = 0$ over `[0, 2π]³`.
- Manufactured solution `cos(x) sin(2y) sin(3z)` integrates to 0.
- Pairwise and Kahan agree within `1e-13` relative on a non-trivial-scale
  integrand.
- Bit-equality across repeated invocations (placeholder for the
  `OMP_NUM_THREADS`-varying test that lands in Phase 20).

### Documentation

- New User Guide note `docs/user-guide/notes/phase-3-domain-integration.md`.
- New Developer Note `docs/developer-note/notes/phase-3-domain-integration.md`,
  first to use the five-section + external-citation rule introduced in
  `specs/CLAUDE.md` step 4d (commits `aabfbf6`, `2309ed4`).

## 0.2.0 — Phase 2: gradient and divergence (scalar)

### Added

- `gridcalc::diff::gradient(const Field<double>&) -> Field<Vec3d>` — 2nd-order
  central-difference gradient on a periodic scalar field. Returns Cartesian
  triple `(∂f/∂x, ∂f/∂y, ∂f/∂z)` per grid point.
- `gridcalc::diff::divergence(const Field<Vec3d>&) -> Field<double>` —
  2nd-order central-difference divergence on a periodic vector field.
- `gridcalc::stencil::FirstDerivative<Order>` — sibling to
  `Coefficients<Order>` (which is the second-derivative table). `Order = 2`
  specialization with weights `{-0.5, 0, 0.5}` (consumer scales by `1/h`).
  Primary template intentionally undefined; `Order = 4` deferred to Phase 7.
- `gridcalc/core/EigenAliases.hpp` — new home for project-local Eigen
  typedefs (`Vec3d` today; `Mat3d` and tensor aliases will land here when
  Phase 13/14 needs them).

### Changed

- `Vec3d` alias relocated from `core/Grid.hpp` to `core/EigenAliases.hpp`.
  Fully-qualified name `gridcalc::core::Vec3d` is unchanged; Phase 1 callers
  recompile without source changes.

### Tests

- Convergence-order tests pass for `gradient` and `divergence` on trig inputs
  over `N ∈ {16, 32, 64}` (log-log slope in `[1.8, 2.2]`, ratio in
  `[3.5, 4.5]`).
- Round-trip identity test: `divergence(gradient(ψ))` converges at order 2
  to the analytical Laplacian `-(kx²+ky²+kz²) ψ` on the same sweep.

## 0.1.0 — Phase 1: periodic 3D scalar grid + Laplacian

### Added

- `gridcalc::core::Grid` — Cartesian-orthogonal sampling grid with per-axis
  cell sizes; i-fastest storage layout.
- `gridcalc::core::Field<T, Policy = Periodic>` — 3D field on a `Grid` with a
  pluggable out-of-range index policy. Three constructors (zero-init,
  broadcast value, callable sampler).
- `gridcalc::core::IndexPolicy::Periodic` — wrap-around index policy.
- `gridcalc::stencil::Coefficients<Order>` (second-derivative table),
  `Order = 2` specialized; `Order = 4` deferred to Phase 7.
- `gridcalc::diff::laplacian` — 2nd-order central-difference periodic
  Laplacian.

### Tests

- Trig-eigenvalue recovery and 2nd-order convergence pass on `N ∈ {16, 32, 64}`.

## 0.0.1 — Phase 0: project scaffold

### Added

- CMake 3.20+ project skeleton, Eigen 3.4.0 + GoogleTest 1.14.0 pinned,
  `.clang-format` / `.clang-tidy` / `CMakePresets.json`, GitHub Actions CI on
  Ubuntu (GCC + Clang). Proprietary `LICENSE` (all rights reserved).
