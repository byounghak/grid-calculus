# Validation: Phase 10 — Documentation infrastructure

## Build

- [ ] CI green on Ubuntu GCC and Ubuntu Clang (test jobs unchanged from Phase 9).
- [ ] New `docs-build` CI job green: produces three artifacts in the `gridcalc-docs-pdfs` bundle — Doxygen output (tarball or directory), `user-guide.pdf`, `developer-note.pdf`.
- [ ] New `docs-lint` CI job green: `scripts/check-since.py` finds no missing `\since` tags; Doxygen with `WARN_AS_ERROR=YES` produces no warnings.
- [ ] Old `Render Doc PDFs` CI job removed; `scripts/render-docs.sh` deleted; `gridcalc-docs-pdfs` artifact bundle continues to exist (now produced by the new job, with the new contents).
- [ ] Local build succeeds: from a clean checkout with TeX Live + Doxygen + Graphviz installed, `cmake --preset clang-debug -DGRIDCALC_BUILD_DOCS=ON && cmake --build build/clang-debug --target gridcalc-docs-doxygen gridcalc-docs-userguide gridcalc-docs-developernote` produces all three artifacts. `scripts/build-docs.sh` does the same.
- [ ] Pandoc removed from the project's required toolchain (was added post-Phase-4 for the lightweight render). `specs/tech-stack.md` updated to match.
- [ ] `gridcalc` CMake target remains INTERFACE.
- [ ] `clang-format` and `clang-tidy` pass on any newly added headers (this phase is mostly docs/build/CI; no new headers expected).

## Tests

- [ ] All Phase 1–9 tests still pass without modification (134 FD tests + 5 spectral basic + 70 FD-FFT cross-check = 209 tests).
- [ ] Version test updated: `VersionTest.ReturnsZeroTenZero` passes against `gridcalc::version() == "0.10.0"`.

## Documentation

- [ ] `docs/Doxyfile` configured for public headers (`include/gridcalc/`, excluding `detail/`); `WARN_IF_UNDOCUMENTED=YES`, `WARN_AS_ERROR=YES`; `HAVE_DOT=YES`; `GENERATE_LATEX=YES`.
- [ ] User Guide LaTeX skeleton (`docs/user-guide/main.tex`, `memoir` class) builds cleanly via `latexmk -pdf` and produces `user-guide.pdf` containing:
  - Front matter: title, copyright, TOC.
  - Chapter 01: "Hello grid" tutorial covering Phases 1–9 in a single end-to-end walkthrough.
  - Chapters 02–10: nine absorbed phase notes (Phase 1 through Phase 9), one chapter per phase.
  - Chapters 11–12: placeholder stubs for Phase 12 (Cahn–Hilliard) and Phase 16 (diamond-lattice).
- [ ] Developer Note LaTeX skeleton (`docs/developer-note/main.tex`, `book` class) builds cleanly via `latexmk -pdf` and produces `developer-note.pdf` containing:
  - Front matter: title, copyright, TOC.
  - Chapter 00: overview / how-to-read-this-book.
  - Chapters 01–09: nine absorbed phase notes, one chapter per phase, each preserving the five-section structure (Theory · Math derivation · Algorithm · Design decisions · References).
  - Back matter: auto-generated Doxygen API reference via `\input{build/docs/doxygen/latex/refman.tex}`. No broken cross-references between chapters and the API reference.
- [ ] `docs/{user-guide,developer-note}/notes/` directories are empty (the per-phase markdown files are migrated). Each directory's `README.md` updated to redirect future content to the LaTeX chapters.
- [ ] `docs/README.md` documents:
  - Local build commands for all three artifacts (cmake-driven + `scripts/build-docs.sh`).
  - Required local toolchain (TeX Live + Doxygen + Graphviz).
  - Fallback to `scripts/get-docs.sh --ci` for fetching CI-built PDFs without a local toolchain.
- [ ] Every public symbol shipped through Phase 9 carries a `\since` tag. Verified by both `scripts/check-since.py` and Doxygen's `WARN_IF_UNDOCUMENTED`.
- [ ] `CHANGELOG.md` entry under `0.10.0` lists the Doxyfile, LaTeX skeletons, chapter migration, retired pandoc render, new docs-build CI job, and `\since` lint.
- [ ] `STATUS.md` updated: Phase 10 row marked Done; Phase 11 stated as the new Next action. New "Decisions worth knowing" entries for the full LaTeX toolchain, the emptied `notes/` directories, and the `\since` lint as a CI gate.
- [ ] `specs/tech-stack.md` updated: pandoc removed from required toolchain; `\since` lint policy noted.

## Performance

- [ ] **Build-cost budget: < 90 s added to a clean CI run.** Measured on the CI run for the merge commit.
  - `docs-build` job (parallel with test jobs): expected ~30–60 s wall-clock for Doxygen + `latexmk` × 2.
  - `docs-lint` job (parallel): expected ~5–10 s.
  - Net new time on the critical path is the longer of the two (< 90 s comfortably).
