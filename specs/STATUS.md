# STATUS

Hand-off snapshot. Update this file whenever a phase completes or a major decision changes.

**Last updated:** 2026-05-01

## Where the project stands

**Phase 1 — Periodic 3D scalar grid + Laplacian is done.** Version bumped to `0.1.0`. The library now has real numerics: `gridcalc::core::Grid` (Cartesian-orthogonal, per-axis cell sizes), `gridcalc::core::Field<T, Policy=Periodic>` (i-fastest storage, three constructors including callable-init), `gridcalc::core::IndexPolicy::Periodic`, `gridcalc::stencil::Coefficients<2>` (with the `Order=4` extension point hooked but unspecialized), and `gridcalc::diff::laplacian` returning a fresh `Field<double>`. Acceptance tests pass: trig eigenvalue recovered to relative max-norm `~9.6e-3` at `N=32`, log-log convergence slope ~2 across `N ∈ {16, 32, 64}`. CMake target `gridcalc` is still INTERFACE (header-only); Eigen is propagated as a SYSTEM include.

**Previously (Phase 0):** PR #1 merged into `main` on 2026-05-01 with the buildable empty skeleton — CMake 3.20+ project, Eigen 3.4.0 + GoogleTest v1.14.0 pinned in `cmake/Dependencies.cmake`, repo layout per `tech-stack.md`, `.clang-format` / `.clang-tidy` / `CMakePresets.json`, GitHub Actions CI on Ubuntu (GCC + Clang), and the proprietary LICENSE.

**Repository:** [github.com/byounghak/grid-calculus](https://github.com/byounghak/grid-calculus) — **private**. Local `main` tracks `origin/main` over SSH (`git@github.com:byounghak/grid-calculus.git`). Two commits pushed:

- `docs: add initial constitution and workflow specs` (`mission`, `tech-stack`, `roadmap`, `CLAUDE`, `STATUS`)
- `chore: add .gitignore for build artifacts, doc outputs, editor clutter`

**License:** Proprietary. All rights reserved by the author; no open-source license. The `LICENSE` file added in Phase 0 will carry an "All rights reserved" notice rather than a permissive license. Mission/scope unchanged: production / industrial use, distributed to authorized recipients only.

## What's in `specs/`

- [mission.md](mission.md) — purpose, audience (production/industrial materials science), quality bar, in/out-of-scope.
- [tech-stack.md](tech-stack.md) — C++17 + Eigen + GoogleTest + Google Benchmark + PocketFFT (verification only) + OpenMP. CMake 3.20+. Docs: Doxygen + LaTeX User Guide (memoir) + Developer Note (book) → PDFs only (no Sphinx, no public HTML site since the repo is private).
- [roadmap.md](roadmap.md) — 23 phases, 0 through 22. Phase 0 = project scaffold; Phase 21 = cross-platform CI hardening; Phase 22 = v1.0 release.
- [CLAUDE.md](CLAUDE.md) — feature workflow runbook. Auto-loaded by Claude Code in this directory.
- `STATUS.md` — this file.

## Next action

**Phase 2 — Gradient and divergence (scalar).** Per `CLAUDE.md`, the entry point is:

1. Open branch `YYYY-MM-DD-phase-2-gradient-divergence` (use today's date when starting).
2. Use `AskUserQuestion` to pin down API choices (`gradient` return type — `Field<Vec3d>`? Three `Field<double>` components? — naming for the divergence input type, whether to lift `Vec3d` out of `core/Grid.hpp` into a shared `core/EigenAliases.hpp` now that Phase 2 needs it).
3. Create `specs/YYYY-MM-DD-phase-2-gradient-divergence/{plan,requirements,validation}.md`.
4. Implement: `diff/Gradient.hpp`, `diff/Divergence.hpp`, plus the round-trip identity test `divergence(gradient(ψ)) ≈ laplacian(ψ)` and convergence-order tests for both operators on trig inputs.

**Acceptance** (from `roadmap.md` Phase 2): convergence-order tests pass for both `gradient` and `divergence` on trig inputs.

A scheduled follow-up agent (`trig_01M1NGo52vJhJEjx4NyvoEdf`) will check in on 2026-05-22 to decide whether to file a tracking issue for Phase 21.

## Phase progress

| Phase  | Title                                                   | Status      |
| ------ | ------------------------------------------------------- | ----------- |
| 0      | Project scaffold                                        | Done        |
| 1      | Periodic 3D scalar grid + Laplacian                     | Done        |
| 2–9    | Remaining periodic FD operators + spectral verification | Not started |
| 10     | Documentation infrastructure                            | Not started |
| 11–14  | Higher-order functionals, CH demo, vector/tensor fields | Not started |
| 15–17  | Lattice basis, multi-atom basis, sublattice operators   | Not started |
| 18–19  | Non-periodic BCs, implicit diffusion                    | Not started |
| 20     | Performance pass                                        | Not started |
| 21     | Cross-platform CI hardening (Apple Clang + MSVC)        | Not started |
| 22     | v1.0 release                                            | Not started |

## Decisions worth knowing (encoded in the specs but easy to miss)

- **Gradients are reported in Cartesian** regardless of lattice basis. Lattice geometry is an implementation detail of `Grid`. (Locked in Phase 15.)
- **Unified mixing-by-default differentiation on multi-atom-basis grids.** At a site of sublattice α, gradients pull contributions from all nearby sites regardless of sublattice. The decoupled per-sublattice mode is opt-in via a policy tag. Weights are derived per `(α, derivative, accuracy)` by generalized finite differences (Taylor-moment matching). (Phase 16.)
- **Functional API for high-rank data** must use contraction expression templates (Phases 11 and 14) to avoid materializing rank-6 intermediate tensors (3⁶ = 729 components/point for a 4th-order gradient of rank-2 data).
- **FFT is verification only** — never the production differentiation engine, never a solver. PocketFFT, header-only, behind `-DGRIDCALC_USE_FFT=ON` (default on). The CI cross-check between FD and FFT operators is a permanent gate from Phase 9 onward.
- **Documentation toolchain split.** User Guide = `memoir` class (typography); Developer Note = `book` class because `doxygen.sty` conflicts with `memoir` on `\hangpara` / `\hangparas` and several preamble packages. `pdflatex` is locked. **Sphinx was dropped** when the repo went private — PDFs are now the only narrative-documentation output, and the "Hello grid" content moved into the User Guide as a chapter (was originally planned as a Sphinx page).
- **Repository visibility: private; license: proprietary.** No open-source license; redistribution requires authorization. `LICENSE` file added in Phase 0 will be an "All rights reserved" notice. There is no public docs hosting; PDFs ship as private-repo release assets only. Mission target ("production / industrial") is unchanged.
- **Phase 0 ships Linux-only CI.** The original roadmap acceptance bar for Phase 0 listed all four compiler families (GCC, Clang, Apple Clang, MSVC). Decided in the Phase 0 spec round to narrow that to Ubuntu GCC + Ubuntu Clang for cheaper iteration, and to insert a new **Phase 21 — Cross-platform CI hardening** before the v1.0 release phase (which moved to Phase 22). The full four-family compiler-support matrix in `tech-stack.md` is unchanged — it just isn't enforced in CI until Phase 21. Known accumulated-portability-bug risk recorded in `roadmap.md` "Risks to watch."
- **`byounghak/atomic-structure` integration deferred to Phase 15+.** Considered importing the existing private library [`atomic-structure`](https://github.com/byounghak/atomic-structure) (v1.6.0; locally at `/Users/byounghak/Projects/KMC refactor/Atomic structure`) as the backing implementation of Phase 1's `Grid`. Decided against: Phase 1 is Cartesian-orthogonal only and a 30-line in-house struct is the right weight; atomic-structure's `Box` + `Crystal` + `Structure` map naturally onto **Phase 15 (non-orthogonal Bravais)** and **Phase 16 (multi-atom basis)** instead. Integration approach (replace `Grid` internals vs. consume `Structure` directly) is decided in the Phase 15 spec round. Recorded as an Implementation note under both Phase 15 and Phase 16 deliverables in `roadmap.md`.
- **Code style.** Leading-underscore member vars (`_lattice`, never `_Capital` or `__double`); camelCase verb-prefixed methods (`getX`, `isY`, `findZ`, `buildW`); snake_case struct fields; PascalCase types; STL-protocol carve-out (`begin`/`end`/`size`/`empty`/`data`/`swap` keep canonical names as 1-line shims that delegate to the verb-prefixed accessor).
- **Compiler tier.** GCC, Clang, Apple Clang, MSVC are all blocking-CI tier-1. The "Windows best-effort" line in `tech-stack.md` is a support-priority statement, not a CI demotion.

## Open / deferred items

- **Eigen unsupported tensor module risk.** May fall back to in-house `Tensor3` / `Tensor4` if `Eigen::TensorFixedSize` doesn't build cleanly across the compiler matrix. Decision deferred until Phase 13 (vector/tensor fields) actually uses it.

## How to resume in a new session

A fresh Claude Code session in this directory auto-loads `CLAUDE.md`. Read this `STATUS.md` first to get current state, then tell Claude: *"Continue from the next phase per STATUS.md."* Claude should run the `CLAUDE.md` runbook from step 1.

When a phase completes, update this file:

1. Move the phase row to "Done" (add a new section if you want to track completion dates).
2. Bump the **Last updated** date.
3. Add any new decisions that emerged to "Decisions worth knowing."
4. Resolve any closed items from "Open / deferred."
