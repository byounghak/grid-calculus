# Plan: Phase 0 — Project scaffold

Four commits, all on branch `2026-05-01-phase-0-project-scaffold`. Order keeps the spec docs landing first so the rest of the work is review-able against an agreed plan, and the global-spec edits land before the implementation that depends on the new acceptance bar.

## Group 1 — Phase 0 spec files

- Add `specs/2026-05-01-phase-0-project-scaffold/plan.md` (this file).
- Add `specs/2026-05-01-phase-0-project-scaffold/requirements.md`.
- Add `specs/2026-05-01-phase-0-project-scaffold/validation.md`.

**Commit message:** `docs: add Phase 0 specs (plan, requirements, validation)`

## Group 2 — Global spec revisions

- `specs/roadmap.md`: revise Phase 0 acceptance to "passes on GCC + Clang on Ubuntu". Insert new **Phase 21 — Cross-platform CI hardening (Apple Clang + MSVC)** before the v1.0 release phase; bump the existing "Phase 21 — Documentation polish, tutorial, v1.0 release" to **Phase 22**. Add a corresponding bullet under "Risks to watch" noting that portability bugs against Apple Clang / MSVC may accumulate until Phase 21.
- `specs/STATUS.md`: bump **Last updated**; rewrite the **Next action** section so it reflects Phase 0 in progress on this branch (rather than "not started"); update the phase-progress table to renumber the row for v1.0 release and insert the new cross-platform phase; add a "Decisions worth knowing" entry for the Linux-only Phase 0 + new Phase 21.
- `specs/tech-stack.md`: in the Repository-layout block, change the `docs/` annotation from "Sphinx sources, Doxygen config" to "Doxygen config + LaTeX sources for User Guide and Developer Note" (Sphinx was dropped per STATUS.md but the line was stale).

**Commit message:** `docs: revise roadmap to defer cross-platform CI to Phase 21`

## Group 3 — Project scaffold (build + skeleton + configs)

- `CMakeLists.txt` at repo root:
  - `cmake_minimum_required(VERSION 3.20)`.
  - `project(gridcalc VERSION 0.0.1 LANGUAGES CXX)`.
  - C++17, no compiler extensions.
  - In-source build guard.
  - Warnings: `-Wall -Wextra -Wpedantic -Wconversion` on GCC/Clang, MSVC equivalents (`/W4` etc.). **Warnings-as-errors only in CI** (controlled by `-DGRIDCALC_WARNINGS_AS_ERRORS=ON`, default `OFF` for local dev, set `ON` in the workflow).
  - Include `cmake/Dependencies.cmake`.
  - Define library target `gridcalc` (INTERFACE for now, since the only public symbol is `version()` exposed via a generated header) — this gets promoted to a real STATIC/SHARED target when Phase 1 adds compiled code.
  - `enable_testing()` and `add_subdirectory(test)`.
- `cmake/Dependencies.cmake`:
  - `FetchContent_Declare(Eigen ...)` pinned to tag `3.4.0` (commit hash recorded in a comment alongside).
  - `FetchContent_Declare(GoogleTest ...)` pinned to tag `v1.14.0` (commit hash recorded).
  - `FetchContent_MakeAvailable(...)` for both.
- `include/gridcalc/version.hpp.in` — `configure_file` template returning `"@gridcalc_VERSION@"` from a `constexpr` `gridcalc::version()`.
- `include/gridcalc/{core,tensor,stencil,diff,spectral,func,solve}/.gitkeep` — empty placeholders.
- `src/.gitkeep`, `bench/.gitkeep`, `examples/.gitkeep`, `docs/.gitkeep`, `cmake/.gitkeep` (the last unnecessary if `cmake/Dependencies.cmake` already lives there).
- `test/CMakeLists.txt` and `test/version_test.cpp`:
  - One `TEST(VersionTest, ReturnsZeroZeroOne)` asserting `gridcalc::version() == std::string_view{"0.0.1"}`.
- `.clang-format` — Google base, line length 100, `IndentWidth: 2`, `PointerAlignment: Left`. A header comment lists project-style overrides.
- `.clang-tidy` — `Checks: 'modernize-*,readability-*,performance-*,bugprone-*,-modernize-use-trailing-return-type'` plus an inline rationale comment.
- `CMakePresets.json` — four configure presets (`gcc-debug`, `gcc-release`, `clang-debug`, `clang-release`) each pointing at the same out-of-source build dir under `build/<preset>`.
- `README.md` stub:
  - One-paragraph description (sourced from `mission.md`).
  - "License: proprietary, all rights reserved. See LICENSE."
  - Build snippet: `cmake --preset gcc-debug && cmake --build --preset gcc-debug && ctest --preset gcc-debug` (and a generic fallback `cmake -B build && cmake --build build && ctest --test-dir build`).
  - Pointer to `specs/` for mission, tech-stack, roadmap, status.
- `LICENSE` — short "All rights reserved" notice per `requirements.md`.

**Commit message:** `build: scaffold Phase 0 project (CMake, repo skeleton, version test)`

## Group 4 — GitHub Actions CI

- `.github/workflows/ci.yml`:
  - Trigger: `pull_request` and `push` to `main`.
  - Two jobs in a matrix: `{ os: ubuntu-latest, compiler: [gcc, clang] }`.
  - Each job:
    - Checks out the repo.
    - Installs the compiler if not preinstalled.
    - Configures with the matching CMake preset, `-DGRIDCALC_WARNINGS_AS_ERRORS=ON`.
    - Builds.
    - Runs `ctest --output-on-failure`.
  - Cache `build/_deps/` so repeat runs skip re-downloading Eigen / GoogleTest.
- README.md — add a CI status badge (pointing at the workflow file in the GitHub repo) once the workflow exists. (Defer if the badge URL needs a workflow run to resolve; OK to land later.)

**Commit message:** `ci: add Linux GCC + Clang GitHub Actions workflow`
