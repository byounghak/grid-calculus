# STATUS

Hand-off snapshot. Update this file whenever a phase completes or a major decision changes.

**Last updated:** 2026-05-03 (Phase 7 implementation; PR open, awaiting CI)

## Where the project stands

**Phase 7 — Higher-order accuracy stencils (4th-order) is done.** Version bumped to `0.7.0`. Adds `gridcalc::stencil::Coefficients<4>` (weights `{-1/12, 4/3, -5/2, 4/3, -1/12}`) and `gridcalc::stencil::FirstDerivative<4>` (weights `{1/12, -2/3, 0, 2/3, -1/12}`); both have radius `2` and truncation error `O(h^4)`. Promotes `diff::laplacian`, `diff::gradient`, and `diff::divergence` to function templates on an `Order` parameter (default `2`); all Phase 1–6 callers (`laplacian(field)`, etc.) keep compiling via the C++17 default-template-argument rule. Phase 7 callers write `laplacian<4>(field)` etc. for the new accuracy. Per-axis-independent unitless weights mean anisotropic per-axis spacing (the kind Phase 1's `Grid` supports) needs no special handling. All seven new tests pass: 4th-order convergence sweep on each operator (`slope ∈ [3.5, 4.5]`), order-2-vs-order-4 absolute-accuracy comparison (order 4 ≥ 10× more accurate at `N = 32`), and weight-table verification. All 53 prior tests still pass without modification.

**Phase 7 PR:** open as of 2026-05-03 — *Phase 7 — Higher-order accuracy stencils*. (PR number assigned at push time; STATUS will be refreshed post-merge.)

**Previously (Phase 6):** PR #9 merged into `main` on 2026-05-03 (rebase merge, commit `2ddfbd0`). Tag-dispatched `solve::integrate(psi, rhs, dt, n_steps, ExplicitEuler{}|RK4{})`; classic 4-stage RK4 implementation with heat-equation CFL `0.6963`; templated `solve::diffuse` reads the per-tag CFL limit. Phase 5's `solve::explicitEuler` preserved as a thin wrapper.

**Previously (Phase 5):** PR #7 merged into `main` on 2026-05-02 (rebase merge, commit `54e8c82`). `solve::explicitEuler` (now a thin wrapper) and `solve::diffuse` with the von-Neumann CFL check.

**Previously (Phase 4):** PR #5 merged into `main` on 2026-05-02 (rebase merge, commit `0e6250f`). `func::evaluate` with SFINAE-detected callable arity for `f(ψ)`, `f(ψ, ∇ψ)`, or `f(ψ, ∇ψ, ∇²ψ)`; eager materialization of only the derivatives the callable actually consumes; Ginzburg–Landau hand-computed test passes within stencil-order error.

**Previously (Phase 3):** PR #4 merged into `main` on 2026-05-02 (rebase merge, commit `b1c136d`). `func::integrate` with `Pairwise` (default) and `Kahan` tag dispatch on a periodic-midpoint quadrature rule. New public namespace `gridcalc::func` and directory `include/gridcalc/func/`. Cross-thread-count determinism is trivially met at single-threaded; Phase 20 expands the placeholder test.

**Previously (Phase 2):** PR #3 merged into `main` on 2026-05-02 (rebase merge, commit `ee93091`). `gradient(Field<double>) -> Field<Vec3d>`, `divergence(Field<Vec3d>) -> Field<double>`, `stencil::FirstDerivative<2>`, `core/EigenAliases.hpp`. Convergence-order tests on trig inputs pass; round-trip `divergence(gradient(ψ))` converges at order 2 to the analytical Laplacian.

**Previously (Phase 1):** PR #2 merged into `main` on 2026-05-01 — `Grid`, `Field<T, Policy=Periodic>`, `IndexPolicy::Periodic`, `stencil::Coefficients<2>`, `diff::laplacian`. Trig eigenvalue recovered to relative max-norm `~9.6e-3` at `N=32`, slope ~2 across `N ∈ {16, 32, 64}`. CMake target `gridcalc` remains INTERFACE (header-only); Eigen is propagated as a SYSTEM include.

**Previously (Phase 0):** PR #1 merged into `main` on 2026-05-01 with the buildable empty skeleton — CMake 3.20+ project, Eigen 3.4.0 + GoogleTest v1.14.0 pinned in `cmake/Dependencies.cmake`, repo layout per `tech-stack.md`, `.clang-format` / `.clang-tidy` / `CMakePresets.json`, GitHub Actions CI on Ubuntu (GCC + Clang), and the proprietary LICENSE.

**Repository:** [github.com/byounghak/grid-calculus](https://github.com/byounghak/grid-calculus) — **private**, SSH remote `git@github.com:byounghak/grid-calculus.git`. Merged PRs to date: **#1** (Phase 0), **#2** (Phase 1), **#3** (Phase 2), **#4** (Phase 3), **#5** (Phase 4), **#6** (PDF render workflow + Phase 1/2 doc-notes backfill), **#7** (Phase 5), **#8** (`get-docs.sh` script), **#9** (Phase 6). Phase 7 PR opens 2026-05-03. Current version `0.7.0`. CI: Ubuntu GCC + Ubuntu Clang on every PR plus a `Render Doc PDFs` job that uploads `gridcalc-docs-pdfs` artifacts (Phase 21 widens to Apple Clang + MSVC).

**License:** Proprietary, all rights reserved (`LICENSE` is the short notice agreed in the Phase 0 spec round). No open-source license; redistribution requires authorization. Mission target unchanged: production / industrial use, distributed to authorized recipients only.

## What's in `specs/`

- [mission.md](mission.md) — purpose, audience (production/industrial materials science), quality bar, in/out-of-scope.
- [tech-stack.md](tech-stack.md) — C++17 + Eigen + GoogleTest + Google Benchmark + PocketFFT (verification only) + OpenMP. CMake 3.20+. Docs: Doxygen + LaTeX User Guide (memoir) + Developer Note (book) → PDFs only (no Sphinx, no public HTML site since the repo is private).
- [roadmap.md](roadmap.md) — 23 phases, 0 through 22. Phase 0 = project scaffold; Phase 21 = cross-platform CI hardening; Phase 22 = v1.0 release.
- [CLAUDE.md](CLAUDE.md) — feature workflow runbook. Auto-loaded by Claude Code in this directory.
- `STATUS.md` — this file.

## Next action

**Phase 8 — Higher-order spatial derivatives (3rd, 4th).** Per `CLAUDE.md`, the entry point is:

1. Open branch `YYYY-MM-DD-phase-8-higher-order-derivatives` (use today's date when starting).
2. Use `AskUserQuestion` to pin down API choices: how to add `\partial^3` and `\partial^4` operators (new function names like `diff::thirdDerivative` / `diff::fourthDerivative` vs a templated `diff::derivative<Rank, Order>`); whether mixed partials (`\partial^2 / \partial x \partial y`, etc.) belong here or in a later phase; whether to introduce a constexpr Fornberg generator now (Phase 7 deferred this; Phase 8 is the natural decision point with multiple new orders/ranks landing).
3. Create `specs/YYYY-MM-DD-phase-8-higher-order-derivatives/{plan,requirements,validation}.md`.
4. Implement per the chosen API; convergence-order tests for each new operator.

**Acceptance** (from `roadmap.md` Phase 8): convergence-order tests pass for the new derivatives at the targeted accuracy order.

A scheduled follow-up agent (`trig_01M1NGo52vJhJEjx4NyvoEdf`) will check in on 2026-05-22 to decide whether to file a tracking issue for Phase 21.

## Phase progress

| Phase  | Title                                                   | Status      |
| ------ | ------------------------------------------------------- | ----------- |
| 0      | Project scaffold                                        | Done        |
| 1      | Periodic 3D scalar grid + Laplacian                     | Done        |
| 2      | Gradient and divergence (scalar)                        | Done        |
| 3      | Domain integration (∫ over the grid)                    | Done        |
| 4      | Functional evaluation, callable API                     | Done        |
| 5      | Explicit Euler diffusion solver                         | Done        |
| 6      | RK4 + generic time integrator interface                 | Done        |
| 7      | Higher-order accuracy stencils (4th-order)              | Done        |
| 8–9    | Remaining periodic FD operators + spectral verification | Not started |
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
- **Lightweight pandoc-based PDF render is in CI (post-Phase-4 interim).** Decided after Phase 4 merged. `scripts/render-docs.sh` + the CI `render-docs` job render the markdown notes to two aggregate PDFs (`user-guide.pdf`, `developer-note.pdf`) on every PR/push, uploaded as workflow artifacts. The toolchain is intentionally narrower than Phase 10's plan (`pandoc` + `pdflatex`, no `memoir`/`book`, no Doxygen integration). Phase 10 retires this script and CI job in favor of the full LaTeX skeleton; until then, the lightweight PDFs are the canonical reader-facing artifact.
- **Implicit schemes (Backward Euler, Crank–Nicolson) stay on Phase 19's schedule.** Considered pulling them forward after Phase 6 merged, weighing two alternatives: (a) FFT-diagonalization solver, periodic-only, ~days of work, leverages PocketFFT pin from `tech-stack.md`; (b) matrix-free Conjugate Gradient with a `LinearOperator` abstraction, periodic-or-not, ~week of work, sets up infrastructure that Phase 19 will need anyway for non-periodic BCs. Decided to defer — keep the original roadmap order, do Phase 7 (4th-order accuracy stencils) next. Phase 19 will introduce implicit schemes from scratch when Phase 18 lands non-periodic boundary policies; option (b) is the natural starting point at that point. Recorded so the trade-off doesn't get re-litigated when Phase 19 starts.
- **Per-phase User Guide / Developer Note notes from Phase 1 forward.** Decided after Phase 2 merged (rule applied from Phase 3 onward at first); extended to cover Phases 1 and 2 retroactively after Phase 4 merged, when the lightweight pandoc PDF render landed and the user wanted Phases 1 and 2 included in the resulting PDFs. Every feature PR from Phase 1 forward has a paired markdown note under `docs/user-guide/notes/phase-N-<slug>.md` (user-facing API summary + worked example) and `docs/developer-note/notes/phase-N-<slug>.md`. The Developer Note follows a fixed five-section structure (**Theory · Math derivation · Algorithm · Design decisions · References**) per `specs/CLAUDE.md` step 4d, with at least one external citation per note (DOI URL, full bibliographic entry, or permanent URL). `validation.md` carries an explicit checkbox enforcing both the structure and the non-empty References section. Phase 10 absorbs the accumulated notes into the LaTeX chapters and empties the two `notes/` directories.

## Open / deferred items

- **Eigen unsupported tensor module risk.** May fall back to in-house `Tensor3` / `Tensor4` if `Eigen::TensorFixedSize` doesn't build cleanly across the compiler matrix. Decision deferred until Phase 13 (vector/tensor fields) actually uses it.
- **Boundary-condition tests for the Euler diffusion solver.** Deferred to Phase 18. Phase 5's `solve::diffuse` and `solve::explicitEuler` are periodic-only by construction (only `IndexPolicy::Periodic` exists); the Phase 18 Dirichlet / Neumann work needs to add the corresponding solver tests (Gaussian decay against an absorbing wall, no-flux mass conservation at a reflective wall, etc.) at the same time as the new boundary policies. Tracked here so it isn't lost when Phase 18 starts.

## How to resume in a new session

A fresh Claude Code session in this directory auto-loads `CLAUDE.md`. Read this `STATUS.md` first to get current state, then tell Claude: *"Continue from the next phase per STATUS.md."* Claude should run the `CLAUDE.md` runbook from step 1.

When a phase completes, update this file:

1. Move the phase row to "Done" (add a new section if you want to track completion dates).
2. Bump the **Last updated** date.
3. Add any new decisions that emerged to "Decisions worth knowing."
4. Resolve any closed items from "Open / deferred."
