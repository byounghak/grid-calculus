# STATUS

Hand-off snapshot. Update this file whenever a phase completes or a major decision changes.

**Last updated:** 2026-05-02 (Phase 2 merged into `main`)

## Where the project stands

**Phase 2 — Gradient and divergence (scalar) is done.** Version bumped to `0.2.0`. The library now exposes the basic first-order operators: `gridcalc::diff::gradient(Field<double>) -> Field<Vec3d>` and `gridcalc::diff::divergence(Field<Vec3d>) -> Field<double>`, both 2nd-order central differences applied per axis with periodic wrap delegated to the input `Field`'s `Policy`. The first-derivative stencil lives in a new sibling template `gridcalc::stencil::FirstDerivative<Order>` (with `Order=2` specialized), separate from the Phase 1 `Coefficients<Order>` second-derivative table. The `Vec3d` alias was lifted out of `core/Grid.hpp` into a new shared `core/EigenAliases.hpp` (fully-qualified name unchanged). All seven new tests pass: trig-gradient recovery, gradient convergence sweep (slope ~2 on `N ∈ {16,32,64}`), analytical-divergence recovery on `V = (sin x cos y, cos x sin y, sin z)`, divergence convergence sweep on the same V, and the round-trip identity `divergence(gradient(ψ))` converging at order 2 to the analytical Laplacian.

**Phase 2 PR:** **#3** — *Phase 2 — Gradient and divergence (scalar)* — merged into `main` on 2026-05-02 (rebase merge, commit `ee93091`). CI was green on Ubuntu GCC and Ubuntu Clang.

**Previously (Phase 1):** PR #2 merged into `main` on 2026-05-01 — `Grid`, `Field<T, Policy=Periodic>`, `IndexPolicy::Periodic`, `stencil::Coefficients<2>`, `diff::laplacian`. Trig eigenvalue recovered to relative max-norm `~9.6e-3` at `N=32`, slope ~2 across `N ∈ {16, 32, 64}`. CMake target `gridcalc` remains INTERFACE (header-only); Eigen is propagated as a SYSTEM include.

**Previously (Phase 0):** PR #1 merged into `main` on 2026-05-01 with the buildable empty skeleton — CMake 3.20+ project, Eigen 3.4.0 + GoogleTest v1.14.0 pinned in `cmake/Dependencies.cmake`, repo layout per `tech-stack.md`, `.clang-format` / `.clang-tidy` / `CMakePresets.json`, GitHub Actions CI on Ubuntu (GCC + Clang), and the proprietary LICENSE.

**Repository:** [github.com/byounghak/grid-calculus](https://github.com/byounghak/grid-calculus) — **private**, SSH remote `git@github.com:byounghak/grid-calculus.git`. Merged PRs to date: **#1** (Phase 0 — project scaffold), **#2** (Phase 1 — periodic 3D scalar grid + Laplacian), **#3** (Phase 2 — gradient and divergence). Current version `0.2.0`. CI: Ubuntu GCC + Ubuntu Clang on every PR (Phase 21 widens to Apple Clang + MSVC).

**License:** Proprietary, all rights reserved (`LICENSE` is the short notice agreed in the Phase 0 spec round). No open-source license; redistribution requires authorization. Mission target unchanged: production / industrial use, distributed to authorized recipients only.

## What's in `specs/`

- [mission.md](mission.md) — purpose, audience (production/industrial materials science), quality bar, in/out-of-scope.
- [tech-stack.md](tech-stack.md) — C++17 + Eigen + GoogleTest + Google Benchmark + PocketFFT (verification only) + OpenMP. CMake 3.20+. Docs: Doxygen + LaTeX User Guide (memoir) + Developer Note (book) → PDFs only (no Sphinx, no public HTML site since the repo is private).
- [roadmap.md](roadmap.md) — 23 phases, 0 through 22. Phase 0 = project scaffold; Phase 21 = cross-platform CI hardening; Phase 22 = v1.0 release.
- [CLAUDE.md](CLAUDE.md) — feature workflow runbook. Auto-loaded by Claude Code in this directory.
- `STATUS.md` — this file.

## Next action

**Phase 3 — Domain integration (∫ over the grid).** Per `CLAUDE.md`, the entry point is:

1. Open branch `YYYY-MM-DD-phase-3-domain-integration` (use today's date when starting).
2. Use `AskUserQuestion` to pin down API choices (reduction strategy — pairwise vs. Kahan, deterministic-across-thread-counts contract, naming under `func/Integrate.hpp`, whether to expose a public `integrate` overload taking a callable for the integrand or only `integrate(Field<double>)`).
3. Create `specs/YYYY-MM-DD-phase-3-domain-integration/{plan,requirements,validation}.md`.
4. Implement: `func/Integrate.hpp` with deterministic reduction order; tests for `∫ 1 dV = domain volume` and `∫ sin(kx) dV = 0` over a period; thread-count-determinism test.

**Acceptance** (from `roadmap.md` Phase 3): reduction is deterministic across thread counts (verified by test).

A scheduled follow-up agent (`trig_01M1NGo52vJhJEjx4NyvoEdf`) will check in on 2026-05-22 to decide whether to file a tracking issue for Phase 21.

## Phase progress

| Phase  | Title                                                   | Status      |
| ------ | ------------------------------------------------------- | ----------- |
| 0      | Project scaffold                                        | Done        |
| 1      | Periodic 3D scalar grid + Laplacian                     | Done        |
| 2      | Gradient and divergence (scalar)                        | Done        |
| 3–9    | Remaining periodic FD operators + spectral verification | Not started |
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
- **`Field<T, Policy = IndexPolicy::Periodic>` is templated on policy from Phase 1.** Only `Periodic` exists today; `Dirichlet` and `Neumann` arrive at Phase 18. The policy template parameter was the *one* piece of forward-looking design admitted at Phase 1 specifically so that Phase 18 lands as an additive change rather than a public-API break. Everything else in Phase 1 stays YAGNI.
- **`Field<T>` storage layout is i-fastest** (`linear_index = i + Nx * (j + Ny * k)`) and is part of the public contract from Phase 1 forward. Matches FFT libraries (PocketFFT, FFTW) and Eigen's column-major default. Phase 9 spectral verification and Phase 20 benchmarks both depend on this; a layout flip later would invalidate Phase-20 baselines.
- **`stencil::Coefficients<Order>` is a template from Phase 1**, with `Order = 2` specialized and `Order = 4` deferred to Phase 7. The primary template is intentionally undefined so unsupported orders trip at compile time at the call site rather than silently zeroing the stencil. Phase 7 adds a specialization with no call-site refactor.
- **First-derivative coefficients live in a sibling template `stencil::FirstDerivative<Order>` (Phase 2), not in `Coefficients<Order>`.** Reusing the `Coefficients` name for first-derivative weights would silently change the semantics of Phase 1's existing call sites (Laplacian uses second-derivative coefficients). Two narrowly-named tables alongside each other was the clean call. A unifying redesign — for example a 2-parameter `Coefficients<Derivative, Order>` — is deferred to Phase 7, when both tables gain `Order = 4` specializations and the trade-offs are easier to evaluate.
- **`core/EigenAliases.hpp` is the home for project-local Eigen typedefs.** Created in Phase 2 when `gradient` and `divergence` joined the Phase 1 `Grid` as call sites for `Vec3d`. Phase 1's `Vec3d` was relocated here without changing its fully-qualified name (`gridcalc::core::Vec3d`). Future tensor/matrix aliases (`Mat3d` etc., as needed by Phase 13/14) belong here as well.
- **`divergence(gradient(ψ))` is *not* identically equal to `laplacian(ψ)` at finite `h`.** The composed operator is effectively a stride-2 second-derivative stencil (`(f(i+2) - 2 f(i) + f(i-2)) / (4 h²)`) whereas `laplacian` uses the 3-point stencil `(f(i-1) - 2 f(i) + f(i+1)) / h²`. Both converge to `∇²` at order 2 but with different leading constants; Phase 2's round-trip test therefore compares to the analytical Laplacian, not to `laplacian()`. This is decided behavior, not a bug.
- **Eigen propagated as a `SYSTEM` include.** `gridcalc`'s `INTERFACE_INCLUDE_DIRECTORIES` for Eigen are added with the `SYSTEM` keyword in the root `CMakeLists.txt`, so `-Wconversion` (and the rest of our strict warning set) stays active on our code but does not fire on `<Eigen/...>`. Without this, `WARNINGS_AS_ERRORS=ON` in CI breaks on Eigen's NEON typecasting headers.
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
