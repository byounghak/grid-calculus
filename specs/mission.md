# Mission

## What this library is

A C++ library for computing **functional values and their derivatives on data defined over a structured crystal-lattice grid**. The primary use cases are phase-field models, continuum field theories, and reaction–diffusion systems on crystalline domains.

Given a field $\psi$ (scalar, vector, or rank-2 tensor) sampled on a structured 3D lattice, the library evaluates expressions of the form

$$
F[\psi] \;=\; \int_\Omega f\bigl(\psi,\;\nabla\psi,\;\nabla^2\psi,\;\ldots,\;\nabla^4\psi\bigr)\,dV
$$

and supports time integration of diffusion-type PDEs

$$
\partial_t \psi \;=\; \mathcal{L}[\psi]
$$

where $\mathcal{L}$ may involve gradients up to fourth order (e.g. Cahn–Hilliard, phase-field crystal).

## Scope

### In scope

- **Structured 3D grids on a Bravais lattice** with arbitrary (orthogonal or non-orthogonal) primitive vectors $\mathbf{a}_1,\mathbf{a}_2,\mathbf{a}_3$ — e.g. simple cubic, FCC, BCC, hexagonal, triclinic.
- **Multi-atom unit cells** (lattice + basis): each Bravais cell may contain $M$ sites at specified fractional positions $\{\boldsymbol{\tau}_\alpha\}_{\alpha=1}^{M}$ — e.g. diamond, zincblende, rocksalt (NaCl), HCP, graphene. Fields may be defined per sublattice.
- **Field data of rank 0, 1, 2** (scalar / vector / tensor) at each grid point.
- **Spatial derivatives up to order 4** via finite-difference stencils with selectable accuracy.
- **Spectral (FFT-based) reference derivatives** on periodic grids — used as the gold-standard truth against which every finite-difference operator is automatically cross-checked in CI. The FD path remains the production code; the FFT path is for verification.
- **Boundary conditions**: periodic first; Dirichlet, Neumann, mixed in later phases.
- **Functional evaluation**: integration of user-supplied densities $f$ over the domain.
- **Diffusion-type time integration**: explicit and semi-implicit schemes.
- **Performance**: cache-friendly memory layout, OpenMP shared-memory parallelism.

### Out of scope (initially)

- Unstructured / adaptive meshes.
- GPU / distributed-memory parallelism (MPI). May be revisited if a clear need emerges.
- Spectral (FFT-based) methods as the **primary** differentiation engine, and FFT-based **diffusion solvers**. The FFT path exists only as a verification harness for the periodic FD operators; the production differentiation and time-stepping code is finite-difference-based.
- Equation systems beyond diffusion class (e.g. wave, hyperbolic, full elasticity).
- Built-in optimization, inverse problems, or AD over the functional.

## Audience and quality bar

Target audience: **production / industrial use** by computational materials scientists and engineers. This sets a high bar:

- **API stability.** Public interfaces follow [Semantic Versioning](https://semver.org/). Breaking changes only at major versions, with migration notes.
- **Correctness.** Every numerical operator has a unit test against an analytic answer (e.g. derivatives of $\sin(kx)$, polynomials) confirming the expected order of accuracy. Convergence tests run in CI.
- **Documentation.** Every public class and free function has Doxygen-style docs. The repository ships with a tutorial that walks from "hello grid" to a running phase-field example.
- **Performance.** Each released version reports benchmark numbers for the canonical operators (Laplacian, biharmonic, functional integration) on a fixed hardware baseline. Regressions block release.
- **Reproducibility.** Builds are deterministic given a pinned compiler and dependency set. Floating-point reduction order is documented where it matters.
- **Long-term support.** No silent removals. Deprecations carry at least one minor-version grace period before removal.

## Non-goals

- Being the fastest possible solver for any single problem; we trade some peak speed for generality and clarity.
- Replacing dedicated phase-field codes (MOOSE, PRISMS-PF). We are a **library** that such codes can be built on top of, not an end-user simulation framework.
- Supporting every conceivable boundary condition or grid topology. We commit to a small set and do them well.

## Success criteria

- A user can express a phase-field free-energy functional $F[\psi]$ in code that reads close to the math, evaluate it on a periodic 3D grid, and time-step its gradient flow — all using documented, stable APIs.
- All shipped operators converge at their advertised order on manufactured solutions.
- Every periodic finite-difference operator is independently cross-checked against an FFT-based spectral reference in CI; the FD-vs-FFT discrepancy scales at the FD's advertised order.
- A single-precision-equivalent runtime for the canonical Laplacian on a 256³ grid is published and tracked release-over-release.
