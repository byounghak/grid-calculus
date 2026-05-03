# Tech Stack

This document pins the technology choices for the library. Changes require a written rationale and a version bump per the policy in `mission.md`.

## Language

- **C++17.** Chosen for broadest compiler reach (including older HPC clusters and vendor toolchains common in industrial settings) while still providing `if constexpr`, structured bindings, `std::optional`, `std::variant`, and parallel algorithms.
- We **do not** use compiler-specific extensions in public headers. Internal `.cpp` files may use `[[gnu::...]]`-style attributes guarded by feature macros.
- Where a C++20/23 feature would be valuable (concepts, `std::mdspan`, `std::span`), we provide a small internal shim and revisit when the language baseline moves.

## Compiler & platform support

- Linux primary target (HPC). macOS supported. Windows best-effort.

## Numerical core

- **[Eigen 3.4+](https://eigen.tuxfamily.org/)** for:
  - Fixed-size small types: `Eigen::Vector3d`, `Eigen::Matrix3d`.
  - Rank-3 and rank-4 small tensors via `Eigen::TensorFixedSize` (from `unsupported/Eigen/CXX11/Tensor`) or a thin in-house wrapper if the unsupported module proves unstable.
  - Sparse matrices and direct/iterative solvers (`SparseLU`, `BiCGSTAB`, `ConjugateGradient`) for implicit time integration on non-periodic grids in later phases.
- **[PocketFFT](https://gitlab.mpcdf.mpg.de/mtr/pocketfft)** for FFT-based spectral reference derivatives on periodic grids. Header-only, BSD-3-Clause-licensed, no transitive system dependencies. **Vendored** under `extern/pocketfft/pocketfft_hdronly.h` rather than `FetchContent`-fetched (decided in the Phase 9 spec round): the project's mission target is "production / industrial," for which build determinism without network access is valuable. The vendored revision is recorded in `extern/pocketfft/VERSION.txt`. The spectral path is **not** the production differentiation engine — it exists as an independent reference against which the periodic finite-difference operators are cross-checked in CI (see `specs/roadmap.md`, Phase 9). Available behind `gridcalc::spectral::*`. Gated by `-DGRIDCALC_USE_FFT=ON/OFF` (default `ON`); turning it off disables the spectral verification tests but does not affect the rest of the library.

## Parallelism

- **OpenMP** for shared-memory parallelism on stencil sweeps and reductions. Required to be optional at build time (`-DGRIDCALC_USE_OPENMP=ON/OFF`) so single-threaded builds remain valid.
- No MPI, no CUDA, no SYCL in the initial scope. Hooks are kept clean (no global state in operators) so a future portability layer is feasible without rewriting kernels.

## Build system

- **CMake 3.20+**. Targets exported as `gridcalc::gridcalc` for `find_package(gridcalc)`. Version, compile features, and include directories are propagated through usage requirements; downstream users do not set flags by hand.
- **Out-of-source builds only.** In-source builds are blocked by a CMake guard.
- **Compiler support matrix** (CI-tested):
  - GCC 9, 11, 13
  - Clang 13, 16
  - Apple Clang 14+
  - MSVC 19.30+ (Visual Studio 2022)
- Warnings: `-Wall -Wextra -Wpedantic -Wconversion` (and MSVC equivalents) treated as errors in CI. Library headers must compile cleanly under these flags for downstream consumers.

## Testing

- **GoogleTest** for unit and integration tests. Chosen over Catch2 for its mature death-test, parameterized-test, and fixture story, which we will lean on for convergence-order tests parameterized by stencil order and grid resolution.
- Every public operator ships with:
  1. A correctness test on an analytic input (polynomial / trigonometric).
  2. A convergence-order test confirming the asymptotic rate (e.g. 2nd-order stencil → slope 2 in log-log).
- **CI** runs the full test matrix on every PR. A nightly job runs longer convergence sweeps and benchmarks.

## Benchmarking

- **[Google Benchmark](https://github.com/google/benchmark)** for microbenchmarks of stencil kernels and reductions. Results are written to a JSON artifact per release; the release script diffs against the previous version and flags regressions over a documented threshold.

## Static analysis and formatting

- **clang-format** enforced in CI. Style file lives at the repo root; deviations from upstream Google/LLVM style are documented inline.
- **clang-tidy** with a curated `.clang-tidy` (modernize, readability, performance, bugprone). Findings are warnings in CI today, errors at v1.0.

## Code style

- Class member variables use a **leading underscore** (e.g., `_lattice`, `_input_box`). Names always start with `_` followed by a lowercase letter — never `_Capital` or `__double_underscore`, both of which are reserved by the C++ standard.
- Class member methods (instance and static) are **camelCase and start with a verb** (e.g., `getLattice`, `getCellCounts`, `getCommensurateBox`, `findIndex`, `isDefined`, `buildByDistance`). Conventions: `get*` for accessors / computed reads, `is*` for boolean predicates, `find*` for lookups that may throw, `build*` / `compute*` / etc. for active operations. `operator()` / other operator overloads are exempt; non-method names (constructors, struct fields, type aliases) follow their own conventions (snake_case fields, PascalCase types).
- **STL-protocol carve-out.** Methods that the standard library looks up by name to enable interop — `begin`, `end`, `cbegin`, `cend`, `rbegin`, `rend`, `size`, `empty`, `data`, `swap` — keep their canonical lowercase STL names so the class participates in range-for, `std::ranges`, structured bindings, `std::swap`, etc. When a class exposes such a method, the verb-prefixed name (e.g., `getSize`) remains the canonical project-style accessor and the STL-named method (e.g., `size`) is added as a 1-line shim that delegates to it.

## Documentation

- **Doxygen** for API reference, generated from header comments. Graphviz `dot` enables class / include / call graphs.
- Header comments use Doxygen `///` triple-slash markers with backslash tags (`\brief`, `\param`, `\returns`, `\throws`, `\note`, `\file`, `\since`). Every public class, function, member variable, and type alias in the public headers carries a Doxygen comment, including a `\since` tag noting the version it was introduced.
- **Developer Note** and **User Guide** are authored in LaTeX and built to PDF (infrastructure stands up in Phase 10; final v1.0 content polish lands in Phase 22). PDFs are the canonical narrative documentation; there is no hosted HTML site (the repository is private and we do not maintain a public docs mirror). The Developer Note is hand-written and `\input{}`s Doxygen-generated LaTeX subfiles for its auto-built API reference and Graphviz diagrams; the User Guide is fully hand-authored with no Doxygen integration.
- **Per-phase doc notes (placeholder, pre-Phase-10).** Until the LaTeX skeletons land in Phase 10, every phase from Phase 3 onward adds two markdown notes — `docs/user-guide/notes/phase-N-<slug>.md` and `docs/developer-note/notes/phase-N-<slug>.md` — in the same PR as the implementation. The User Guide note captures user-facing API additions with a worked example. The Developer Note follows a fixed five-section structure (**Theory · Math derivation · Algorithm · Design decisions · References**) per `specs/CLAUDE.md` step 4d: theory and explicit math derivations sit alongside the algorithm description, and the References section requires at least one external citation per note (DOI URL, full bibliographic entry, or permanent URL). Phase 10's deliverable list includes absorbing the accumulated notes into the corresponding `.tex` chapters and emptying the `notes/` directories.
- **Lightweight pre-Phase-10 PDF render (interim).** Introduced after Phase 4. The script `scripts/render-docs.sh` runs `pandoc --pdf-engine=xelatex` against the accumulated markdown notes and produces two aggregate PDFs — `build/docs/user-guide.pdf` and `build/docs/developer-note.pdf` — each with a TOC, numbered sections, and clickable links. The CI job `render-docs` runs this on every PR/push and uploads the PDFs as artifacts (downloadable from the run page; not committed). This is intentionally a **simpler toolchain** than Phase 10's plan — `pandoc` over plain `xelatex`, no `memoir` or `book` class, no Doxygen integration — so the PDFs are usable now without pulling Phase 10 forward. The engine is `xelatex` (not `pdflatex`) because the markdown notes contain Greek letters and other Unicode glyphs both inside and outside math mode; xelatex handles Unicode natively while pdflatex needs explicit inputenc setup. Phase 10's full LaTeX skeleton stays on `pdflatex` — its hand-authored `.tex` sources can quote Greek via LaTeX commands and benefit from `pdflatex`'s tighter ecosystem with `doxygen.sty`. Phase 10 replaces this lightweight render with the full LaTeX skeleton; until then, the lightweight PDFs are the canonical reader-facing artifact.
- LaTeX toolchain: **`pdflatex`** for both documents (locked early in the documentation work). Document classes are split: the User Guide uses **`memoir`** for nicer reader-facing typography; the Developer Note uses **`book`** because `doxygen.sty` (which the Developer Note `\input{}`s for per-class API subfiles) was designed for `book`-class output and conflicts with memoir on `\hangpara`, `\hangparas`, and several preamble packages (`geometry`, `fancyhdr`, `tocloft`, `changepage`). The split was discovered during smoke testing of the integrated build and is the minimum-cost workaround — both documents share the engine, only the class differs. Driven by `latexmk -pdf`; CI installs `texlive-latex-recommended texlive-latex-extra texlive-fonts-recommended texlive-pictures latexmk` (Ubuntu names) plus `doxygen graphviz`. An engine switch (xelatex / lualatex) requires a documented exception, not a re-vote.
- PDFs are CI build artifacts; `.tex` sources, figures, `.bib`, and the `Doxyfile` live in the repo, but `*.pdf` and Doxygen output directories are gitignored.
- `README.md` and `specs/` for prose; no separate docs hosting site initially.

## Dependency policy

- All dependencies are **header-only or vendor-able**. Transitive system-library requirements beyond a C++17 toolchain are forbidden.
- Dependencies are pulled via CMake's `FetchContent` with **pinned commit hashes**, not version ranges. Updates are explicit changes reviewed in a PR.
- Current pinned set:
  - Eigen (tag `3.4.0` or later patch)
  - GoogleTest (tag `v1.14.0` or later)
  - Google Benchmark (tag `v1.8.x`)
  - PocketFFT — **vendored** under `extern/pocketfft/`, not `FetchContent`. Revision recorded in `extern/pocketfft/VERSION.txt`. Decided in the Phase 9 spec round (build determinism without network access).
- A `cmake/Dependencies.cmake` file is the single source of truth for versions.

## Repository layout

```
gridcalc/
├── CMakeLists.txt
├── cmake/                    Toolchain helpers, Dependencies.cmake, install rules
├── include/gridcalc/         Public headers (installed)
│   ├── core/                 Grid, Field, IndexPolicy, BoundaryCondition
│   ├── tensor/               Small tensor types and contractions
│   ├── stencil/              Stencil definitions and application
│   ├── diff/                 gradient, divergence, laplacian, biharmonic, mixed partials
│   ├── spectral/             FFT-based reference derivatives (verification path, periodic only)
│   ├── func/                 Functional evaluation and integration
│   └── solve/                Time integrators and diffusion solvers
├── src/                      Non-header implementations (kept small)
├── test/                     GoogleTest unit and convergence tests
├── bench/                    Google Benchmark microbenchmarks
├── examples/                 Self-contained programs (Cahn–Hilliard, etc.)
├── docs/                     Doxygen config + LaTeX sources for User Guide and Developer Note
└── specs/                    mission.md, tech-stack.md, roadmap.md (this directory)
```

## Versioning and release

- **Semantic Versioning 2.0.0.** `MAJOR.MINOR.PATCH`.
- Pre-1.0 (`0.y.z`): minor bumps may break API; patch bumps may not.
- Post-1.0: breaking changes only on major; features on minor; bugfixes on patch.
- A `CHANGELOG.md` is updated with every PR that touches public API or behavior.

## What is explicitly **not** in the stack

- No Boost (avoids a heavy transitive dependency that complicates industrial deployment).
- No Conan / vcpkg as required tooling. The CMake `FetchContent` path must always work standalone.
- No Python in the runtime path. A future Python binding (pybind11) is a separate package, not a build-time requirement of this library.
