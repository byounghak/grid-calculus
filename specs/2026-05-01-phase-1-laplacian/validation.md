# Validation: Phase 1 — Periodic 3D scalar grid + Laplacian

The PR does not merge with any unticked items. "Acceptance" rows are pulled from the roadmap Phase 1 bullets; "standard" rows cover hygiene that applies to every phase.

## Build

- [ ] CI green on Ubuntu GCC and Ubuntu Clang (Phase-0 acceptance bar; full four-family bar is still Phase 21).
- [ ] Build is warning-clean under `-Wall -Wextra -Wpedantic -Wconversion` with `GRIDCALC_WARNINGS_AS_ERRORS=ON`.
- [ ] `clang-format -n --Werror` clean on every changed `.hpp` / `.cpp` (run locally; CI enforcement still deferred).

## Acceptance (from roadmap)

- [ ] `LaplacianTest.RecoversTrigEigenvalue` passes — Laplacian of $\sin(k_x x)\sin(k_y y)\sin(k_z z)$ at $\mathbf{k} = (1, 2, 1)$, $N=32$, $L=2\pi$ recovers $-6$ × the input within relative max-norm `5e-3`.
- [ ] `LaplacianTest.SecondOrderConvergence` passes — sweep $N \in \{16, 32, 64\}$ on $L=2\pi$, $\mathbf{k}=(1,1,1)$. Error factor between successive resolutions in `[3.5, 4.5]`; least-squares log-log slope in `[1.8, 2.2]`.

## Tests (new)

- [ ] `test/index_policy_test.cpp` covers in-range, negative, and large-overflow `Periodic::wrap`.
- [ ] `test/grid_test.cpp` covers shape exposure, cell volume, linear-index formula, and rejection of zero/negative N or non-positive cell sizes.
- [ ] `test/field_test.cpp` covers all three constructors (zero, value, callable), the linear-index round-trip via `operator()`, and periodic wrap on negative and overflow indices.
- [ ] `test/laplacian_test.cpp` includes the two acceptance tests above plus the output-grid-shape sanity check.
- [ ] `ctest --output-on-failure` exits zero locally on at least one of `clang-debug` / `gcc-debug` and in CI on both.

## Documentation

- [ ] Every new public symbol in `core/Grid.hpp`, `core/Field.hpp`, `core/IndexPolicy.hpp`, `stencil/CentralDifference.hpp`, `diff/Laplacian.hpp` carries a Doxygen block (`\brief`, `\param`, `\returns` where applicable, `\since 0.1.0`).
- [ ] No narrative tutorial added — Phase 10 is responsible for docs infrastructure; Phase-1 content folds into the User Guide's "Hello grid" chapter at Phase 10. Decision recorded here so a future reviewer doesn't flag the omission.
- [ ] `CHANGELOG.md` — **still not required at Phase 1**; first CHANGELOG entry lands when there is a public API decision worth communicating beyond the version bump (likely Phase 2 or 4). Decision recorded here.

## Performance

- [ ] **Not applicable at Phase 1.** No benchmarks; no `bench/baselines/` entries.

## Workflow hygiene

- [ ] All five planned commits land in order with the prefixes from `plan.md`.
- [ ] STATUS.md and CMakeLists version bump land in the **same** commit (Group 5) so a reader checking out Phase-1 mid-merge never sees a stale STATUS or version.
- [ ] No `Co-Authored-By: Claude` trailer on any commit.
- [ ] PR opened against `main`; CI watched until both jobs report green; merge only after every checkbox above is ticked.
