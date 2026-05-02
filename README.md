# gridcalc

Differential operators and functional evaluation on periodic crystalline grids,
in C++17. Aimed at production / industrial materials-science workflows; see
[`specs/mission.md`](specs/mission.md) for the full statement of purpose.

## Status

Pre-alpha. **Phase 0 — Project scaffold** is the current work; see
[`specs/STATUS.md`](specs/STATUS.md) for what has landed and what is next.

## License

**Proprietary.** All rights reserved. Distribution is to authorized recipients
only. See [`LICENSE`](LICENSE) for the full notice.

## Build

Requires CMake 3.20+ and a C++17 compiler. Eigen and GoogleTest are fetched
automatically via CMake `FetchContent` on first configure (no system install
required).

Using CMake presets (recommended):

```sh
cmake --preset gcc-debug
cmake --build --preset gcc-debug
ctest --preset gcc-debug
```

Or without presets:

```sh
cmake -B build && cmake --build build && ctest --test-dir build
```

In-source builds are blocked by a CMake guard.

## Documentation

The authoritative project documents live under [`specs/`](specs/):

- [`mission.md`](specs/mission.md) — purpose, audience, quality bar, scope.
- [`tech-stack.md`](specs/tech-stack.md) — language, dependencies, build, code style.
- [`roadmap.md`](specs/roadmap.md) — 23-phase plan from scaffold to v1.0.
- [`STATUS.md`](specs/STATUS.md) — current state, decisions worth knowing, next action.
- [`CLAUDE.md`](specs/CLAUDE.md) — feature workflow runbook.

Per-feature spec directories under `specs/YYYY-MM-DD-feature-name/` carry the
plan, requirements, and validation checklist for each phase.

API reference (Doxygen) and the User Guide / Developer Note PDFs are stood up
in Phase 10 of the roadmap; they are not yet built in this scaffold.
