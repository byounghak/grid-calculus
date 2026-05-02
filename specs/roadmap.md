# Roadmap

Phases are intentionally **small** — each is a 1–3 day chunk landing as a single reviewable PR with tests. Every phase ends in a green CI build and a tagged pre-release.

The ordering is deliberate: we drive a thin slice (scalar field, periodic, basic ∇ and ∇²) to **functional evaluation** as quickly as possible to validate the integration story, then widen scope.

## Legend

- **Goal** — what changes about the library after this phase.
- **Deliverables** — concrete artifacts (code + tests + docs).
- **Acceptance** — the bar that must be met to call the phase done.

---

## Phase 0 — Project scaffold

- **Goal.** A buildable, testable empty library.
- **Deliverables.**
  - `CMakeLists.txt`, `cmake/Dependencies.cmake` pinning Eigen + GoogleTest.
  - Repository layout per `tech-stack.md`.
  - One trivial test (`gridcalc::version()` returns a known string).
  - GitHub Actions CI: Ubuntu GCC + Ubuntu Clang, both green. (Apple Clang and MSVC jobs are deferred to Phase 21 — Cross-platform CI hardening.)
  - clang-format and clang-tidy configs.
  - `README.md` stub linking to `specs/`.
- **Acceptance.** Fresh clone → `cmake -B build && cmake --build build && ctest --test-dir build` passes on Ubuntu GCC and Ubuntu Clang. The full four-family CI bar from `tech-stack.md`'s compiler-support matrix is enforced at Phase 21, not here.

## Phase 1 — Periodic 3D scalar grid + Laplacian

- **Goal.** First real numerics: a scalar field on a periodic 3D grid with a 2nd-order-accurate Laplacian.
- **Deliverables.**
  - `core/Grid.hpp` — Nx/Ny/Nz, lattice vectors, cell volume. Cartesian-orthogonal only for now.
  - `core/Field.hpp` — `Field<double>` with linear indexing.
  - `core/IndexPolicy.hpp` — `Periodic` policy with wrapping `index(i,j,k)`.
  - `stencil/CentralDifference.hpp` — 2nd-order central differences in 1D.
  - `diff/Laplacian.hpp` — `laplacian(field) -> Field<double>`.
- **Acceptance.**
  - Apply Laplacian to $\sin(k_x x)\sin(k_y y)\sin(k_z z)$ → recover $-(k_x^2+k_y^2+k_z^2)$ times the input.
  - Convergence test: halving $h$ reduces error by 4×, log-log slope ≈ 2.

## Phase 2 — Gradient and divergence (scalar)

- **Goal.** Complete the basic first-order operators.
- **Deliverables.**
  - `diff/Gradient.hpp` — `gradient(field) -> Field<Vector3d>`.
  - `diff/Divergence.hpp` — `divergence(vfield) -> Field<double>`.
  - Round-trip identity test: `divergence(gradient(ψ)) ≈ laplacian(ψ)` to expected order.
- **Acceptance.** Convergence-order test passes for both operators on trig inputs.

## Phase 3 — Domain integration (∫ over the grid)

- **Goal.** A primitive that turns a scalar field into a single number — the integration backbone of functional evaluation.
- **Deliverables.**
  - `func/Integrate.hpp` — `integrate(Field<double>) -> double` with deterministic reduction order (pairwise or Kahan).
  - Test: `∫ 1 dV` returns domain volume; `∫ sin(kx) dV` returns 0 over a period.
- **Acceptance.** Reduction is deterministic across thread counts (verified by test).

## Phase 4 — Functional evaluation, callable API (scalar, periodic)

- **Goal.** First end-to-end functional: $F[\psi] = \int f(\psi, \nabla\psi, \nabla^2\psi)\,dV$ with user-supplied $f$.
- **Deliverables.**
  - `func/Functional.hpp` — accepts a callable `f(double, Vector3d, double)` and evaluates pointwise then integrates.
  - Worked example: Ginzburg–Landau free energy $F = \int \tfrac{1}{2}|\nabla\psi|^2 + W(\psi^2-1)^2 \,dV$.
  - Tutorial doc page walking through the example.
- **Acceptance.** Numerical $F$ matches a hand-computed reference for a chosen $\psi$ field within stencil-order error.

## Phase 5 — Explicit Euler diffusion solver

- **Goal.** Time integration enters the picture: $\partial_t\psi = D\nabla^2\psi$.
- **Deliverables.**
  - `solve/ExplicitEuler.hpp` — fixed time-step integrator parameterized by an RHS callable.
  - `solve/Diffusion.hpp` — convenience driver `diffuse(psi, D, dt, n_steps)`.
  - Stability check: warns/errors when $D\,dt/h^2 > 0.5/d$ (CFL for $d$-D explicit diffusion).
- **Acceptance.** Numerical solution to a Gaussian initial condition matches the analytic solution within finite-difference error after 100 steps.

## Phase 6 — RK4 + generic time integrator interface

- **Goal.** Cleaner API for time-stepping; higher-order time accuracy.
- **Deliverables.**
  - `solve/Integrator.hpp` interface and `RK4` implementation.
  - Refactor `Diffusion` to consume the integrator interface.
- **Acceptance.** Convergence-in-time test: error scales as $\Delta t^4$.

## Phase 7 — Higher-order accuracy stencils (4th-order)

- **Goal.** The accuracy-switching mechanism, validated on existing operators.
- **Deliverables.**
  - 4th-order central-difference coefficients for 1st and 2nd derivatives.
  - `Order` template / runtime parameter on `gradient`, `laplacian`.
  - Convergence tests for both orders.
- **Acceptance.** Log-log convergence slopes match advertised orders within tolerance.

## Phase 8 — Higher-order spatial derivatives (3rd, 4th)

- **Goal.** Mixed partials and the biharmonic operator.
- **Deliverables.**
  - `diff/MixedPartial.hpp` — $\partial^2/\partial x_i\partial x_j$.
  - `diff/Biharmonic.hpp` — $\nabla^4$.
  - 4th-order accuracy variants of both.
- **Acceptance.** Convergence tests pass; biharmonic of $\sin(\mathbf{k}\cdot\mathbf{x})$ recovers $|\mathbf{k}|^4$.

## Phase 9 — Spectral verification harness

- **Goal.** A second, independent path for computing derivatives on **periodic** grids — used as the gold-standard reference against which every finite-difference operator is automatically validated. The FFT path is **not** exposed as a primary differentiation engine; the FD operators remain the production code.
- **Scope at this phase.** Orthogonal Cartesian periodic grids only. Non-orthogonal-lattice and multi-atom-basis spectral counterparts are added as deliverables of Phases 15 and 16 alongside the corresponding grid features.
- **Deliverables.**
  - `spectral/Fft.hpp` — thin wrapper over PocketFFT for forward/backward 3D real-to-complex FFTs on `Field<double>`.
  - `spectral/Derivatives.hpp` — spectral counterparts of the periodic FD operators delivered through Phase 8:
    - `spectral::gradient(field) -> Field<Vector3d>`
    - `spectral::laplacian(field) -> Field<double>`
    - `spectral::biharmonic(field) -> Field<double>`
    - `spectral::mixed_partial(field, i, j) -> Field<double>`

    Implemented as multiplication by $i\mathbf{k}$, $-|\mathbf{k}|^2$, $|\mathbf{k}|^4$, etc. in Fourier space, then inverse transform.
  - **Wavenumber and aliasing discipline.** Standard convention $k_x = 2\pi n/(N_x h_x)$ for $n \in [-N/2, N/2)$. The Nyquist mode is zeroed for odd-order derivatives to avoid aliasing artifacts. Documented in the operator headers.
  - **CI cross-check fixture.** A parameterized GoogleTest fixture that, for every periodic FD operator at every advertised accuracy order, applies the FFT version to the same input and verifies the FD–FFT discrepancy scales as $\mathcal{O}(h^p)$ where $p$ is the FD order. Becomes a permanent CI gate; future FD changes that silently break accuracy are caught here.
  - **Build flag.** `-DGRIDCALC_USE_FFT=ON` (default). Disabling skips the spectral tests but does not affect FD code.
  - `bench/spectral_vs_fd.cpp` — reports per-operator runtime ratio for documentation.
- **Acceptance.**
  - **Spectral correctness.** Spectral $\nabla^2 \sin(\mathbf{k}\cdot\mathbf{x}) = -|\mathbf{k}|^2 \sin(\mathbf{k}\cdot\mathbf{x})$ to relative error $< 10^{-12}$ on a $64^3$ grid.
  - **Round-trip identity.** `ifft(fft(ψ)) ≈ ψ` to within $10\,\epsilon_\text{mach}$ on a random $64^3$ field.
  - **FD-vs-FFT cross-check.** For each of `gradient`, `laplacian`, `biharmonic`, `mixed_partial` at 2nd- and 4th-order accuracy, the FD–FFT difference scales at the FD's advertised order across at least three resolutions. Verified in CI.
  - **Build-cost budget.** Enabling `GRIDCALC_USE_FFT` adds < 5 s to a clean build (PocketFFT is small and header-only).

## Phase 10 — Documentation infrastructure

- **Goal.** Stand up the full documentation toolchain so every later phase can land its own incremental docs (API reference, User Guide chapters, Developer Note sections) into a working build. This phase is plumbing — each subsequent feature phase remains responsible for the doc deliverables it already lists (Cahn–Hilliard tutorial in Phase 12, diamond-lattice example in Phase 16, etc.); they fold in here rather than waiting for the v1.0 polish at Phase 22.
- **Deliverables.**
  - `docs/Doxyfile` configured for the public headers; Graphviz `dot` enabled for class / include / call graphs.
  - LaTeX skeletons for the **Developer Note** (`book` class; `\input{}`s Doxygen-generated subfiles for the API reference and Graphviz diagrams) and the **User Guide** (`memoir` class; fully hand-authored). Both built with `latexmk -pdf`.
  - A "Hello grid" tutorial chapter in the User Guide covering the operators delivered through Phase 9 (scalar periodic FD, integration, callable functional, explicit/RK4 diffusion, higher-order accuracy, biharmonic + mixed partials, FFT verification).
  - `docs/README.md` documenting the local build commands for each artifact.
  - CI job: installs `texlive-latex-recommended texlive-latex-extra texlive-fonts-recommended texlive-pictures latexmk doxygen graphviz`; produces three artifacts on every PR — Doxygen output, Developer Note PDF, User Guide PDF. PDFs are not committed.
  - **`\since` policy enforcement.** A CI lint scans the public headers and fails the build if any public class, function, member variable, or type alias is missing a `\since` tag. All symbols delivered through Phase 9 carry tags before this phase merges.
- **Acceptance.**
  - All three documentation artifacts build cleanly in CI on every PR.
  - The `\since` lint blocks PRs that introduce public symbols without the tag.
  - Build-cost budget: full docs build adds < 90 s to a clean CI run.
  - The User Guide PDF contains the "Hello grid" chapter; the Developer Note PDF contains an auto-built API reference for everything in Phases 0–9 with no broken cross-references.

## Phase 11 — Functional with up-to-4th-order gradients

- **Goal.** Extend the functional API to densities depending on $\nabla^3$ and $\nabla^4$ terms.
- **Deliverables.**
  - Updated `func/Functional.hpp` accepting callables with extended argument lists.
  - Worked example: phase-field-crystal free energy $F = \int \tfrac{1}{2}\psi[(q_0^2+\nabla^2)^2 - \epsilon]\psi + \tfrac{1}{4}\psi^4 \,dV$.
- **Acceptance.** PFC functional value matches a hand-computed reference.

## Phase 12 — Cahn–Hilliard example end-to-end

- **Goal.** A real, scientifically meaningful demo proving the layered abstractions hold up.
- **Deliverables.**
  - `examples/cahn_hilliard.cpp` solving $\partial_t\psi = M\nabla^2(\delta F/\delta\psi)$.
  - VTK output of $\psi$ at user-specified intervals.
  - Tutorial page: equation → code → results, with screenshots.
- **Acceptance.** Coarsening dynamics qualitatively match published CH benchmarks.

## Phase 13 — Vector and tensor fields

- **Goal.** Generalize `Field<T>` beyond scalars; rank-raising gradients.
- **Deliverables.**
  - `Field<Vector3d>`, `Field<Matrix3d>` instantiations.
  - `gradient(Field<Vector3d>) -> Field<Matrix3d>` (∂_i v_j).
  - `divergence(Field<Matrix3d>) -> Field<Vector3d>`.
  - `tensor/Contraction.hpp` — single and double contractions.
- **Acceptance.** Vector identities (e.g. $\nabla\times\nabla\phi = 0$) recovered to discretization error.

## Phase 14 — Functional evaluation for vector/tensor data

- **Goal.** Functionals over vector and tensor fields.
- **Deliverables.**
  - Generalized `Functional` that accepts callables over `(T, Gradient<T>, ...)` argument packs.
  - Worked example: linear elastic energy $F = \int \tfrac{1}{2}\,\boldsymbol{\sigma}:\boldsymbol{\varepsilon}\,dV$.
- **Acceptance.** Elastic energy of a simple stretched cubic block matches the analytical value.

## Phase 15 — Non-orthogonal Bravais lattices

- **Goal.** Support general single-atom Bravais lattices (FCC, BCC, hexagonal, triclinic) by introducing a basis transform between Cartesian and lattice coordinates.
- **Deliverables.**
  - `Grid` carries primitive vectors $\mathbf{a}_1,\mathbf{a}_2,\mathbf{a}_3$ and a Cartesian-↔-lattice transform $A = [\mathbf{a}_1\ \mathbf{a}_2\ \mathbf{a}_3]$.
  - Stencils continue to operate on lattice indices; derivatives are converted to Cartesian via $A^{-T}$ on the gradient (chain rule).
  - Decision documented: **derivatives are reported in Cartesian coordinates** regardless of grid basis. The user works in physically meaningful units; the lattice geometry is an implementation detail of the grid.
  - Tests: derivative of $\sin(\mathbf{k}\cdot\mathbf{x})$ on a sheared / hexagonal grid matches Cartesian truth.
  - **Implementation note (decided 2026-05-01).** This phase is the planned integration point for [`byounghak/atomic-structure`](https://github.com/byounghak/atomic-structure) (currently v1.6.0 in `/Users/byounghak/Projects/KMC refactor/Atomic structure`). Its `Box` (origin + integer-direction columns + Cartesian lengths), `Crystal` (lattice matrix + basis atoms), and `Structure` (constructed lattice with `getCellCounts()` / `findIndex(basis_i, n1, n2, n3)` / minimum-image distance) map directly onto what this phase needs. The Phase 15 spec round will decide between (a) replacing `Grid`'s internals with `NSPC_Atomic_Structure::Structure` while keeping the public `Grid` API stable, vs. (b) consuming `Structure` directly. The dep enters via `FetchContent` against a tagged release; gridcalc's pinned Eigen must take precedence over atomic-structure's `FindOrFetchEigen.cmake`. Phase 1's `Grid` is intentionally an in-house Cartesian-only struct so that integration happens here rather than at the foundation.
- **Acceptance.** Convergence tests pass on at least one non-orthogonal Bravais lattice (hexagonal recommended as the smoke test).

## Phase 16 — Multi-atom basis (lattice with basis), unified differentiation

- **Goal.** Generalize the data model from "one site per Bravais cell" to "$M$ sites per Bravais cell at fractional offsets $\{\boldsymbol{\tau}_\alpha\}$." Differential operators treat the multi-basis crystal as a **single decorated point lattice**: by default, the derivative at a site mixes contributions from **all** nearby sites regardless of sublattice. Decoupled per-sublattice differentiation is available as an opt-in.
- **Deliverables.**
  - `core/Basis.hpp` — a `Basis` value type holding $M$ and the fractional offsets $\{\boldsymbol{\tau}_\alpha\}$. Validated: offsets in $[0,1)^3$, no duplicates within a tolerance. **Implementation note:** if Phase 15 adopted `byounghak/atomic-structure`, `Basis` should be a thin adapter over `NSPC_Atomic_Structure::BasisAtom` (drop the species/name fields; keep just `fractional_position`) so Phase 15 and Phase 16 share one underlying lattice representation.
  - `Grid` gains a `Basis` member. The single-atom path is preserved as `Basis::trivial()` (the $M=1$ default); earlier phases keep working unchanged.
  - `Field<T>` carries an explicit **sublattice dimension**: storage shape becomes `(Nx, Ny, Nz, M)`. Indexing helpers: `field(i,j,k,alpha)` and `field.sublattice(alpha) -> a Bravais-shaped view`.
  - **Unified geometric stencils.** `core/Stencil.hpp` is extended so a stencil is defined by a set of **Cartesian offset vectors**, not lattice-index offsets. At stencil construction, each offset is resolved to a `(Bravais offset, target α)` pair against the chosen `Grid`+`Basis`. Standard 1D central-difference coefficients no longer apply (the resolved points aren't collinear and aren't equally spaced), so weights are derived by **generalized finite differences**: solve a small linear system that matches the Taylor moments of the requested derivative up to the prescribed accuracy order.
  - **Per-sublattice weight tables.** Because the local geometric neighborhood depends on the source sublattice α (e.g. the nearest neighbors of a basis-1 site in diamond are all basis-2, and vice versa), weights are precomputed and cached per `(α, derivative, accuracy_order)` at `Grid` construction. The construction routine performs a consistency check: for each precomputed weight set, verify that the moment conditions match the requested derivative to working precision; fail loudly otherwise.
  - **Default differentiation behavior.** `gradient(field)`, `divergence(field)`, `laplacian(field)`, `biharmonic(field)`, … use the unified stencil. At site α, contributions come from all sublattices β reached by the stencil. The behavior reduces **exactly** to standard central differences when $M=1$.
  - **Opt-in decoupled mode.** A policy tag, e.g. `gradient(field, decoupled_sublattice)`, restricts the stencil to same-α neighbors only. Useful when sublattices carry physically independent order parameters and coupling enters the model elsewhere.
  - `func/Integrate.hpp` updated: domain integration sums over all $M$ sublattices with per-sublattice volume weights (Bravais cell volume divided by $M$ by default; user-supplied weights for non-uniform partitions).
  - Tutorial: a worked example on the **diamond lattice** (FCC + 2-atom basis at $\boldsymbol{\tau}_1=(0,0,0)$, $\boldsymbol{\tau}_2=(\tfrac{1}{4},\tfrac{1}{4},\tfrac{1}{4})$). Show that the default unified $\nabla\psi$ at a basis-1 site uses its four nearest basis-2 neighbors (and farther shells as needed for accuracy), recovers the Cartesian truth, and degrades gracefully if the user opts into decoupled mode on an asymmetric field.
- **Acceptance.**
  - **Unified-default diamond test.** $\psi(\mathbf{x}) = \sin(\mathbf{k}\cdot\mathbf{x})$ sampled at every basis site (single smooth field across both sublattices) → unified $\nabla\psi$ recovers $\mathbf{k}\cos(\mathbf{k}\cdot\mathbf{x})$ at the prescribed order. Verified at 2nd and 4th order accuracy.
  - **Decoupled-mode test.** Two independent sinusoidal fields on the two sublattices → decoupled gradient produces results that depend only on each sublattice's own field; unified gradient produces the (correctly larger) coupled result.
  - **Single-atom regression.** $M=1$ produces bit-identical numerics and bit-identical performance to Phase 15 on the same problem (verified by hash-comparing output and by a benchmark guard).
  - **Storage parity.** $M=1$ storage overhead is zero relative to Phase 15 (verified by `sizeof` and a benchmark check).
  - **Weight-table sanity.** A unit test fails if any precomputed weight set violates its Taylor moment conditions beyond a documented tolerance.

## Phase 17 — Sublattice-aware non-differential operators

- **Goal.** A small operator family for quantities that are **explicitly between sublattices** but are *not* differentiation. These are distinct from the unified differential operators in Phase 16: they are algebraic / geometric reductions over the basis index, useful for tight-binding-like sums, sublattice-resolved order parameters, and inter-sublattice diffusion couplings.
- **Deliverables.**
  - `diff/BondStencil.hpp` — geometric stencils defined by an explicit list of bond vectors (e.g. the four NN bonds in diamond). Distinct from the differential stencils because no derivative is computed; the stencil simply gathers values along bonds.
  - Operators that work on `Field<T>` of shape `(Nx, Ny, Nz, M)` and return a same-shape field:
    - `nearest_neighbor_sum(field, bonds)` — Σ over given bonds.
    - `sublattice_average(field)` — average over the basis index α.
    - `sublattice_difference(field, alpha, beta)` — pointwise $\psi_\alpha - \psi_\beta$.
  - Tests: on a uniform field, NN-sum returns coordination-number × value; on a known antisymmetric mode, `sublattice_difference` returns the analytic amplitude.
- **Acceptance.** Tested on diamond and on a graphene-style 3D-extruded geometry.

## Phase 18 — Dirichlet and Neumann boundary conditions

- **Goal.** First non-periodic support. Stencils degrade gracefully near boundaries.
- **Deliverables.**
  - `IndexPolicy` variants: `Dirichlet`, `Neumann`, `Mixed`.
  - One-sided / biased stencils for cells within (stencil-width − 1) of a non-periodic boundary.
  - BC-aware `Laplacian` and `Gradient`.
- **Acceptance.** Method-of-manufactured-solutions test on a Dirichlet box matches expected order to within boundary degradation.

## Phase 19 — Implicit / semi-implicit diffusion

- **Goal.** Larger stable timesteps for stiff problems.
- **Deliverables.**
  - `solve/CrankNicolson.hpp` building a sparse system via Eigen's `SparseMatrix`.
  - Solver dispatch: CG for SPD systems, BiCGSTAB otherwise.
- **Acceptance.** Stable for $D\,dt/h^2 \gg 1$ on a benchmark; same accuracy as explicit at matching step sizes.

## Phase 20 — Performance pass and benchmarks

- **Goal.** Pin down the v1 performance baseline.
- **Deliverables.**
  - OpenMP parallelization of stencil kernels and reductions.
  - Google Benchmark suite covering Laplacian, biharmonic, integration, full RK4 step.
  - Published baseline numbers in `bench/baselines/`.
- **Acceptance.** Documented scaling to all CI cores; reproducible numbers committed to the repo.

## Phase 21 — Cross-platform CI hardening (Apple Clang + MSVC)

- **Goal.** Bring the GitHub Actions matrix up to the full four-family bar promised in `tech-stack.md` ("Compiler & platform support" + "Compiler support matrix"). Phase 0 deliberately shipped Linux-only CI to keep early iteration cheap; this phase closes the gap before v1.0.
- **Deliverables.**
  - GitHub Actions CI matrix expanded to include macOS (Apple Clang) and Windows (MSVC 19.30+, Visual Studio 2022). Both run the full configure + build + `ctest` pipeline of every prior phase's tests.
  - Whatever portability fixes are required to make the Apple Clang and MSVC jobs green: header `#include` cleanups, `__attribute__((...))` → `[[gnu::...]]`-with-MSVC-fallback macro shims, MSVC-specific warning suppressions, path-separator and case-sensitivity fixes, etc.
  - Optionally widen GCC and Clang version coverage on the Linux jobs (GCC 9/11/13, Clang 13/16) per the support matrix in `tech-stack.md`. May land in this phase or be deferred to a follow-up depending on the size of the portability fixes.
  - `README.md` updated: build snippet covers macOS and Windows in addition to Linux.
- **Acceptance.**
  - All four compiler families (GCC, Clang, Apple Clang, MSVC) are blocking-CI tier-1 with green builds on the latest tagged release of every prior phase's test suite.
  - No regressions on Linux: GCC and Clang jobs remain green throughout the portability work.
  - `tech-stack.md`'s compiler support matrix is now actually enforced rather than aspirational.

## Phase 22 — Documentation polish, tutorial, v1.0 release

- **Goal.** Release-ready per the production/industrial bar in `mission.md`. Most documentation content has accumulated incrementally since Phase 10; this phase is the final pass over completeness, polish, and the release artifacts. Distribution is private — release assets go to authorized recipients, not a public download page.
- **Deliverables.**
  - **User Guide PDF** (memoir class) complete with worked examples spanning every shipped capability: scalar diffusion, Cahn–Hilliard, phase-field-crystal, vector/tensor-field functionals, multi-atom basis (diamond), non-orthogonal Bravais lattices, and Dirichlet/Neumann boundaries.
  - **Developer Note PDF** (book class) complete, including the auto-built API reference subfiles from Doxygen and the architecture/decisions narrative.
  - Both PDFs attached as release assets on the v1.0.0 tag in the private repository.
  - Migration guide template in place for future major bumps.
  - `CHANGELOG.md` finalized for 1.0.0.
  - All public symbols carry `\since` tags (the lint introduced in Phase 10 has been blocking offenders all along; this phase verifies completeness one more time).
- **Acceptance.** Tagged `v1.0.0`; both PDFs attached as private-repo release assets; CI green; benchmarks run on the tag commit.

---

## What we are deliberately **not** doing in v1

These appear on the radar but are out of scope until post-1.0:

- GPU and MPI back-ends.
- Adaptive mesh refinement.
- FFT-based spectral solvers as a primary path.
- Automatic differentiation through the functional.
- Python bindings (separate package).

## Risks to watch

- **Cross-platform portability drift.** Phase 0 ships Linux-only CI (Ubuntu GCC + Clang) and the four-family bar from `tech-stack.md` is not enforced until Phase 21. Apple Clang / MSVC portability bugs may accumulate across Phases 1–20 and surface in bulk at Phase 21. Mitigation: developers with macOS or Windows machines should spot-check builds locally on milestone phases (4, 9, 16); the Phase 21 fix-up budget should be planned as larger than a typical phase.
- **Eigen's unsupported tensor module.** If `Eigen::TensorFixedSize` proves unstable for our rank-3/4 needs, we fall back to a small in-house `Tensor3` / `Tensor4` (deferred to phase 13+). Detection criterion: any compiler in the support matrix fails to build the tensor module headers cleanly.
- **High-rank intermediate tensors.** A 4th-order gradient of rank-2 data is rank-6 (729 components/point). Functionals must be expressible as contractions that **never materialize** the full tensor — this constrains the API design at phase 11 and 14. Contraction expression templates may be needed.
- **Boundary stencils at high accuracy.** A 4th-order biharmonic needs ~9 points; one-sided stencils near boundaries can be ill-conditioned. Phase 18 may require restricting accuracy in the boundary layer; this is acceptable and will be documented.
- **Generalized-FD weight derivation.** With the unified mixing-by-default differentiation in Phase 16, stencil points are not collinear and not equally spaced, so closed-form 1D central-difference coefficients do not apply. Weights are derived per `(α, derivative, accuracy)` by solving a small Taylor-moment linear system at `Grid` construction. Risks: (1) the moment system can be ill-conditioned for poorly-chosen stencil shells — guard with a condition-number check and fall back to a larger shell; (2) weight tables grow with $M$ and accuracy order — measure cache impact in Phase 20 and shard if needed.
- **Cost of unified differentiation.** Per gradient evaluation, the unified default does $\mathcal{O}(K \cdot M)$ work per site versus $\mathcal{O}(K)$ for the decoupled mode. For phase-field workloads this is acceptable and physically correct; it does shift the perf baseline at Phase 20, which should report numbers separately for $M=1$, $M=2$ (diamond), and $M=4$.
- **Spectral verification on multi-basis grids.** When the FFT verification harness is extended in Phases 15–16 to non-orthogonal lattices and multi-atom bases, the spectral path needs a sublattice-aware structure factor (sum over $\alpha$ of $e^{i\mathbf{k}\cdot\boldsymbol{\tau}_\alpha}$) so that derivatives match the FD convention. If this turns out to be subtle or expensive, the cross-check in those configurations may be reduced to a per-sublattice check on smooth fields rather than a full FD-vs-FFT identity test.
