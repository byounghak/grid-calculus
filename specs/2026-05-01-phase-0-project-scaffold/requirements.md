# Requirements: Phase 0 — Project scaffold

## Roadmap reference

> **Goal.** A buildable, testable empty library.
>
> **Deliverables.**
>
> - `CMakeLists.txt`, `cmake/Dependencies.cmake` pinning Eigen + GoogleTest.
> - Repository layout per `tech-stack.md`.
> - One trivial test (`gridcalc::version()` returns a known string).
> - GitHub Actions CI: GCC + Clang + Apple Clang + MSVC, all green.
> - clang-format and clang-tidy configs.
> - `README.md` stub linking to `specs/`.
>
> **Acceptance.** Fresh clone → `cmake -B build && cmake --build build && ctest --test-dir build` passes on all four CI compilers.

Refer to [mission.md](../mission.md) (production/industrial audience, quality bar, scope) and [tech-stack.md](../tech-stack.md) (C++17, Eigen+GoogleTest, code style, repo layout, dependency policy) for global context — this document does not restate them.

## Decisions made for this feature

- **License.** `LICENSE` ships a short "All rights reserved" notice (~5 lines): copyright Byounghak Lee 2026, no use/copy/modify/distribute without prior written permission, "AS IS" warranty disclaimer. Matches STATUS.md's proprietary-distribution stance; not an OSS license.
- **Dependency pin set.** `cmake/Dependencies.cmake` pins **only Eigen and GoogleTest** in Phase 0 (matches the roadmap deliverable text verbatim). Google Benchmark and PocketFFT are pinned in the phases that first consume them (Benchmark → Phase 20, PocketFFT → Phase 9). The "single source of truth" file exists from Phase 0; later phases append entries.
- **CI matrix narrowed to Linux GCC + Clang at Phase 0.** Two jobs: Ubuntu latest with GCC (latest in the apt repos) and Clang (latest). This **does not** satisfy the original four-family acceptance bar in the roadmap, so the bar is being formally revised in this same branch — see `roadmap.md` Phase 0 (revised) and the new Phase 21 (cross-platform CI hardening). STATUS.md is updated to reflect the change.
- **Roadmap structural change.** A new **Phase 21 — Cross-platform CI hardening (Apple Clang + MSVC)** is inserted before the v1.0 release phase, which moves to **Phase 22**. Rationale: the user wants to defer macOS/Windows CI cost until later, but cross-platform support is still a v1.0 prerequisite per `tech-stack.md`'s compiler-support matrix and `mission.md`'s production/industrial bar.
- **Repository layout breadth at Phase 0.** Create all 7 subdirectories under `include/gridcalc/` (`core`, `tensor`, `stencil`, `diff`, `spectral`, `func`, `solve`) tracked via `.gitkeep` markers, plus empty `src/`, `test/`, `bench/`, `examples/`, `docs/`, `cmake/`. Mirrors `tech-stack.md`'s layout from day one so contributors see the intended structure on a fresh clone.
- **Initial version.** `gridcalc::version()` returns `"0.0.1"`. The string is sourced from CMake's `project(... VERSION 0.0.1)` declaration via `configure_file` on `include/gridcalc/version.hpp.in`. Pre-1.0 patch level is appropriate for an empty scaffold per `tech-stack.md`'s "0.y.z: minor bumps may break API" policy.
- **CMakePresets.json.** Ship a small `CMakePresets.json` with `gcc-debug`, `gcc-release`, `clang-debug`, `clang-release` configure presets so dev and CI invoke CMake the same way. Keeps later phases from inventing ad-hoc invocation strings.
- **clang-format style.** Google base style, line length 100, `DerivePointerAlignment: false` + `PointerAlignment: Left`. `IndentWidth: 2`. Code-style overrides documented inline at the top of `.clang-format`.
- **clang-tidy.** `.clang-tidy` enables `modernize-*`, `readability-*`, `performance-*`, `bugprone-*`. Disabled in CI in Phase 0 (added as warnings-only check in a later phase, per `tech-stack.md`'s "warnings in CI today, errors at v1.0" policy). The config file exists so devs can run it locally.
- **`tech-stack.md` hygiene fix in same branch.** Line referencing `docs/` as "Sphinx sources, Doxygen config" is updated to drop Sphinx (already noted in `STATUS.md`'s decisions list but stale in tech-stack.md itself).

## Non-goals

- **No OpenMP wiring at Phase 0.** The `-DGRIDCALC_USE_OPENMP=ON/OFF` option is added in the first phase that needs threading, not now. An empty scaffold has no kernels to parallelize.
- **No PocketFFT or `-DGRIDCALC_USE_FFT` option at Phase 0.** Lands in Phase 9 alongside the spectral verification harness.
- **No Google Benchmark wiring at Phase 0.** Lands in Phase 20.
- **No Doxyfile or LaTeX setup at Phase 0.** All documentation infrastructure is Phase 10's responsibility (per the roadmap). The `docs/` directory exists with a `.gitkeep` only.
- **No install rules at Phase 0.** `install(TARGETS ...)`, `find_package(gridcalc)` config files, and the `gridcalc::gridcalc` exported namespace are wired when there is a real library to install (Phase 1+ has a real library; Phase 0 has only a `version()` shim that doesn't need to be exported yet). This is a small deferral; the roadmap deliverable text doesn't mandate install-tree support at Phase 0.
- **No clang-tidy in CI at Phase 0.** Config file exists; CI invocation deferred.
- **No release tagging at Phase 0.** The roadmap's "every phase ends in a green CI build and a tagged pre-release" line is interpreted as aspirational; explicit pre-release tagging cadence is a project-management decision separate from this feature.

## Deferred questions

- **Apple Clang and MSVC CI jobs** — moved to the new **Phase 21**. Until that phase, portability bugs against those toolchains can accumulate; this is a known and accepted risk noted in `roadmap.md` "Risks to watch."
- **clang-tidy CI enforcement (warnings-only)** — punted to a later phase (likely Phase 1 or whenever the first real kernel lands and there is something for tidy to lint).
- **Install tree, `find_package(gridcalc)` support, and the `gridcalc::gridcalc` exported target** — deferred until there is a real library to consume (Phase 1+).
- **Pre-release tagging cadence** — not pinned in either the roadmap or this spec; will revisit when Phase 1 ships.
