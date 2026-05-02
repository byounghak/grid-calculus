# Validation: Phase 0 — Project scaffold

The PR does not merge with any unticked items. "Acceptance" rows are pulled directly from the (revised) `roadmap.md` Phase 0 deliverables and acceptance bar; "standard" rows cover build hygiene that applies to every phase.

## Build

- [ ] CI green on Ubuntu GCC and Ubuntu Clang (the two blocking jobs at Phase 0; full four-family bar is deferred to Phase 21).
- [ ] CMake configure succeeds with `cmake --preset gcc-debug` and `cmake --preset clang-debug` from a clean checkout.
- [ ] CMake errors out if the user attempts an in-source build.
- [ ] Build is warning-clean under `-Wall -Wextra -Wpedantic -Wconversion` on GCC and Clang (`-DGRIDCALC_WARNINGS_AS_ERRORS=ON` set in CI).
- [ ] `clang-format -n --Werror $(git ls-files '*.hpp' '*.cpp' '*.h' '*.in')` reports clean (run locally; CI enforcement deferred).
- [ ] `.clang-tidy` parses without error (run `clang-tidy --explain-config` locally; CI enforcement deferred).

## Acceptance (from roadmap, revised)

- [ ] Fresh clone → `cmake -B build && cmake --build build && ctest --test-dir build` passes on Ubuntu GCC and Ubuntu Clang. (Fallback "no presets" path explicitly exercised in addition to the preset-based path.)
- [ ] `cmake/Dependencies.cmake` pins Eigen 3.4.0 and GoogleTest v1.14.0 (or later patch versions, per `tech-stack.md` policy) by tag, with the resolved commit hash recorded in a comment alongside.
- [ ] Repository layout matches `tech-stack.md`: `CMakeLists.txt`, `cmake/`, `include/gridcalc/{core,tensor,stencil,diff,spectral,func,solve}/`, `src/`, `test/`, `bench/`, `examples/`, `docs/`, `specs/` all present.
- [ ] `gridcalc::version()` is callable from a downstream `#include <gridcalc/version.hpp>`, returns `"0.0.1"`, and the version string is sourced from CMake's `project(...)` declaration (verified by hand — change `VERSION 0.0.1` to a different value, reconfigure, observe the test fail).
- [ ] `test/` contains exactly one passing `TEST` named `VersionTest.ReturnsZeroZeroOne`.
- [ ] `.clang-format` and `.clang-tidy` both exist at the repo root, both parse cleanly under their respective tools.
- [ ] `README.md` exists, links to `specs/mission.md`, `specs/tech-stack.md`, `specs/roadmap.md`, `specs/STATUS.md`; states "private repository, proprietary license"; includes a build snippet using the CMake presets.
- [ ] `LICENSE` exists with the "All rights reserved" notice per `requirements.md`.

## Tests

- [ ] `VersionTest.ReturnsZeroZeroOne` passes locally on macOS Apple Clang (developer machine spot-check, not CI).
- [ ] `VersionTest.ReturnsZeroZeroOne` passes in CI on both Ubuntu GCC and Ubuntu Clang.
- [ ] `ctest --output-on-failure` exits zero on a clean build.

## Documentation

- [ ] `README.md` build snippet is copy-pastable and works on a fresh clone (verified by running it during local validation).
- [ ] `specs/STATUS.md` "Last updated" bumped to today's date and the Phase-progress table reflects Phase 0 in progress (or done, when this PR merges) and the new Phase 21 row inserted.
- [ ] `specs/roadmap.md` Phase 0 acceptance bar reflects the Linux-only CI scope.
- [ ] `specs/roadmap.md` contains the new Phase 21 — Cross-platform CI hardening section, and the v1.0 release phase is renumbered to Phase 22.
- [ ] `specs/roadmap.md` "Risks to watch" includes the portability-bug-accumulation note.
- [ ] `specs/tech-stack.md` `docs/` annotation no longer mentions Sphinx.
- [ ] `CHANGELOG.md` — **not required at Phase 0** (no public API exists yet); explicit decision recorded here so a future reviewer doesn't flag the omission.

## Performance

- [ ] **Not applicable at Phase 0.** No kernels exist; no benchmarks run; no `bench/baselines/` entries.

## Workflow hygiene

- [ ] All four planned commits land in order (specs → roadmap revisions → scaffold → CI), each with the commit message prefix from `plan.md`.
- [ ] No `Co-Authored-By: Claude` trailer in commit messages (per project-wide preference).
- [ ] No files added that should be `.gitignore`-d (build artifacts, IDE state, OS scratch). Spot-check `git status` before each commit.
