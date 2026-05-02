# STATUS

Hand-off snapshot. Update this file whenever a phase completes or a major decision changes.

**Last updated:** 2026-05-01

## Where the project stands

Pre-implementation. Only the constitution and workflow docs exist in `specs/`; there is no code, no `CMakeLists.txt`, no local git repository yet.

**Repository:** [github.com/byounghak/grid-calculus](https://github.com/byounghak/grid-calculus) — **private**. Remote registered; not yet pushed.

**License:** Proprietary. All rights reserved by the author; no open-source license. The `LICENSE` file added in Phase 0 will carry an "All rights reserved" notice rather than a permissive license. Mission/scope unchanged: production / industrial use, distributed to authorized recipients only.

## What's in `specs/`

- [mission.md](mission.md) — purpose, audience (production/industrial materials science), quality bar, in/out-of-scope.
- [tech-stack.md](tech-stack.md) — C++17 + Eigen + GoogleTest + Google Benchmark + PocketFFT (verification only) + OpenMP. CMake 3.20+. Docs: Doxygen + LaTeX User Guide (memoir) + Developer Note (book) → PDFs only (no Sphinx, no public HTML site since the repo is private).
- [roadmap.md](roadmap.md) — 22 phases, 0 through 21. Phase 0 = project scaffold; Phase 21 = v1.0 release.
- [CLAUDE.md](CLAUDE.md) — feature workflow runbook. Auto-loaded by Claude Code in this directory.
- `STATUS.md` — this file.

## Next action

**Phase 0 — Project scaffold.**

Per `CLAUDE.md`, the entry point is:

1. Open branch `YYYY-MM-DD-phase-0-project-scaffold`.
2. Use `AskUserQuestion` to pin down scaffold-specific details (e.g., GitHub org / repo name, license file, initial namespace name `gridcalc::`).
3. Create `specs/YYYY-MM-DD-phase-0-project-scaffold/` with `plan.md`, `requirements.md`, `validation.md`.
4. Then implement the scaffold: `CMakeLists.txt`, `cmake/Dependencies.cmake` (pinned Eigen + GoogleTest), repo skeleton per `tech-stack.md`, GitHub Actions CI, clang-format / clang-tidy configs, README stub.

**Acceptance** (from roadmap): fresh clone → `cmake -B build && cmake --build build && ctest --test-dir build` passes on GCC + Clang + Apple Clang + MSVC.

## Phase progress

| Phase | Title | Status |
|---|---|---|
| 0 | Project scaffold | Not started |
| 1–9 | Core periodic FD operators + spectral verification | Not started |
| 10 | Documentation infrastructure | Not started |
| 11–14 | Higher-order functionals, CH demo, vector/tensor fields | Not started |
| 15–17 | Lattice basis, multi-atom basis, sublattice operators | Not started |
| 18–19 | Non-periodic BCs, implicit diffusion | Not started |
| 20 | Performance pass | Not started |
| 21 | v1.0 release | Not started |

## Decisions worth knowing (encoded in the specs but easy to miss)

- **Gradients are reported in Cartesian** regardless of lattice basis. Lattice geometry is an implementation detail of `Grid`. (Locked in Phase 15.)
- **Unified mixing-by-default differentiation on multi-atom-basis grids.** At a site of sublattice α, gradients pull contributions from all nearby sites regardless of sublattice. The decoupled per-sublattice mode is opt-in via a policy tag. Weights are derived per `(α, derivative, accuracy)` by generalized finite differences (Taylor-moment matching). (Phase 16.)
- **Functional API for high-rank data** must use contraction expression templates (Phases 11 and 14) to avoid materializing rank-6 intermediate tensors (3⁶ = 729 components/point for a 4th-order gradient of rank-2 data).
- **FFT is verification only** — never the production differentiation engine, never a solver. PocketFFT, header-only, behind `-DGRIDCALC_USE_FFT=ON` (default on). The CI cross-check between FD and FFT operators is a permanent gate from Phase 9 onward.
- **Documentation toolchain split.** User Guide = `memoir` class (typography); Developer Note = `book` class because `doxygen.sty` conflicts with `memoir` on `\hangpara` / `\hangparas` and several preamble packages. `pdflatex` is locked. **Sphinx was dropped** when the repo went private — PDFs are now the only narrative-documentation output, and the "Hello grid" content moved into the User Guide as a chapter (was originally planned as a Sphinx page).
- **Repository visibility: private; license: proprietary.** No open-source license; redistribution requires authorization. `LICENSE` file added in Phase 0 will be an "All rights reserved" notice. There is no public docs hosting; PDFs ship as private-repo release assets only. Mission target ("production / industrial") is unchanged.
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
