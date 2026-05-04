# Requirements: Phase 12 — Cahn–Hilliard example end-to-end

## Roadmap reference

> **Goal.** A real, scientifically meaningful demo proving the layered abstractions hold up.
>
> **Deliverables.**
> - `examples/cahn_hilliard.cpp` solving $\partial_t\psi = M\nabla^2(\delta F/\delta\psi)$.
> - VTK output of $\psi$ at user-specified intervals.
> - Tutorial page: equation → code → results, with screenshots.
>
> **Acceptance.** Coarsening dynamics qualitatively match published CH benchmarks.

(`specs/roadmap.md`, Phase 12.)

## Decisions made for this feature

- **Free-energy form: standard Landau double-well + gradient term.**
  $F[\psi] = \int \!\bigl[-\tfrac12\psi^2 + \tfrac14\psi^4 + \tfrac{\kappa}{2}|\nabla\psi|^2\bigr]\,dV$ with $\kappa = 1$, $M = 1$ in non-dimensional units. The variational derivative is supplied **in closed form** as $\delta F/\delta\psi = -\psi + \psi^3 - \kappa\nabla^2\psi$, hand-coded — no automatic differentiation. The full RHS is therefore $\partial_t\psi = M\nabla^2(-\psi + \psi^3 - \kappa\nabla^2\psi)$, which exercises Phase 1's `diff::laplacian` twice per evaluation and nothing else from later phases of the diff stack.

- **Time integrator: explicit RK4 via Phase 6's `solve::integrate(psi, rhs, dt, n_steps, RK4{})`.** Zero new solver code. The CFL is brutal ($\Delta t \sim h^4/\kappa$ from the $-M\kappa\nabla^4\psi$ linear term) but acceptable on the 32³–64³ demo grid. Semi-implicit / IMEX schemes are explicitly **not** in scope; the project's "FFT is verification only, never a solver" stance (`mission.md`) rules out the spectral-diagonalization path. Implicit treatment of the $\nabla^4$ term remains earmarked for **Phase 19**.

- **VTK output: legacy ASCII `.vtk` STRUCTURED_POINTS, one file per snapshot.** Filename pattern `<prefix>_<step:06d>.vtk`. Hand-written writer (~30–40 lines) under `examples/common/vtk_writer.hpp`; no new build dependency. ParaView and VisIt open the format directly. The writer's interface accepts a `Field<double>`, a `Grid`, an output path, and a field name; it emits `# vtk DataFile Version 3.0`, the `STRUCTURED_POINTS` dataset block (`DIMENSIONS`, `ORIGIN`, `SPACING`), and a `POINT_DATA` block with one `SCALARS` array.

- **Output cadence: per-step interval, runtime CLI flag.** The example's `main` accepts `--snapshot-every N` (default `100`) and writes one VTK file every `N` integration steps. Wall-clock cadence was rejected as adding clock-reading machinery for no scientific gain on a deterministic demo. CLI also exposes `--n-x`, `--dt`, `--n-steps`, `--kappa`, `--mobility`, `--seed`, `--out-dir` (default `./out`).

- **Initial condition: uniform random `±0.05` perturbation around `ψ̄ = 0`,** seeded from a fixed `--seed` (default `42`) so the demo is bit-reproducible. Mass-conservative: the perturbations are zero-mean by construction (we sample `ψ ~ U(-0.05, 0.05)` and then subtract the empirical mean to enforce `Σψ = 0` at machine precision).

- **Box: cubic `L = N·Δx` with `Δx = 1`, periodic on all three axes.** Default `N = 64`; the acceptance test runs at `N = 32` for budget. The interface width `√κ = 1` lets the demo fit ~6–8 emerging domains across each axis at `N = 64`, well into the coarsening regime within reasonable wall-clock.

- **Acceptance bar: visual + monotone-coarsening programmatic check.** Two evidence streams.
  - *Visual.* Three rendered 2D-slice PNGs (early / mid / late times) committed under `docs/user-guide/figures/cahn-hilliard/` and embedded in the User Guide chapter, demonstrating the textbook spinodal-decomposition → coarsening sequence.
  - *Programmatic, in CI.* `test/cahn_hilliard_test.cpp` runs a small CH simulation (`N = 32`, ~3000 RK4 steps, target wall-clock <30 s on Ubuntu Clang) and asserts (a) Lyapunov decay — total free energy `F(t)` is monotonically non-increasing across snapshots within a tolerance set by the explicit-RK4 local truncation error; (b) monotone coarsening — gradient energy $\int|\nabla\psi|^2\,dV$ decreases across post-transient snapshots; (c) mass conservation — $|\bar\psi(t) - \bar\psi(0)| < 10^{-10}$ across the run.
  
  No `R(t) \sim t^{1/3}` power-law fit. The lighter monotone trend is robust against noise at the affordable grid size, and the visual snapshots carry the qualitative-similarity statement.

- **Examples build wiring: new `GRIDCALC_BUILD_EXAMPLES` option, default `OFF`,** mirroring Phase 10's `GRIDCALC_BUILD_DOCS` pattern. Root `CMakeLists.txt` gates `add_subdirectory(examples)` on the option. CI flips the option `ON` on both Ubuntu jobs so the example is type-checked on every PR but not run.

- **RHS extracted to a small reusable header.** `examples/common/cahn_hilliard.hpp` exposes `Field<double> computeCahnHilliardRhs(const Field<double>& psi, double M, double kappa)`. Both the example's `main` and the acceptance test `#include` it; this avoids duplicating the RHS body across the two translation units. Helpers under `examples/` are not part of the library's public API surface and therefore carry **no `\since` tag** — the Phase 10 lint regex covers `include/gridcalc/` only.

- **Snapshot PNG generation.** `scripts/render_ch_snapshots.py` reads the example's `.vtk` output and writes 2D z-midplane slices as PNGs via `matplotlib`. Run locally during this phase; the resulting PNGs are committed. The script is documented in the User Guide chapter so future contributors can regenerate. Not invoked from CI.

## Non-goals

- **No automatic differentiation.** $\delta F/\delta\psi$ is supplied as hand-coded C++. AD over functionals is out of scope for the entire roadmap.
- **No FFT-diagonalized solver.** "FFT is verification only" (`mission.md`) — using the spectral path as a CH integrator violates the locked-in decision and is not on offer.
- **No 2D codepath.** The example is 3D-only (consistent with the library's 3D-first stance from Phase 1). 2D slicing is for visualization, performed in the Python helper post-hoc.
- **No `t^{1/3}` power-law fit in CI.** A separate scaling benchmark could be added later; this phase ships the lighter monotone-coarsening check only.
- **No new public symbols in `include/gridcalc/`.** Phase 12 consumes Phases 1–8 and adds zero library API. The `\since` lint therefore has nothing new to verify; the gate stays green automatically.
- **No new VTK library dependency.** Hand-written writer.
- **No CI invocation of the example.** CI builds the binary (with `GRIDCALC_BUILD_EXAMPLES=ON`) but does not run it; the acceptance test in `test/` provides the runtime gate.

## Deferred questions

- **Implicit / IMEX time-stepping for CH.** Phase 19 (implicit diffusion) — when the matrix-free CG infrastructure lands for non-periodic BCs, an implicit treatment of the $-M\kappa\nabla^4\psi$ term comes naturally with it. Spec round will revisit whether to retrofit `examples/cahn_hilliard.cpp` with an `--integrator imex` switch.
- **`R(t) \sim t^{1/3}` quantitative scaling benchmark.** Could ship as `examples/cahn_hilliard_scaling.cpp` later; not gated, runs longer, fits the log-log slope of `1/k*(t)` or `1/A(t)`. Marked as a follow-up if/when a benchmarking phase wants more rigorous evidence.
- **Vector-valued / multi-component CH.** Phase 13 (vector and tensor fields) opens this up — model B for binary alloys, model H for fluid–fluid mixtures, etc. None in this phase.
- **OpenMP / parallel performance.** Phase 20 (performance pass) — the demo is single-threaded.
