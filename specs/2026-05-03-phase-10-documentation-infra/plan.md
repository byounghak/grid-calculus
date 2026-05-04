# Plan: Phase 10 — Documentation infrastructure

## Group 1 — Doxygen configuration

- Create `docs/Doxyfile` configured for public headers only:
  - `INPUT = include/gridcalc/`.
  - `RECURSIVE = YES`.
  - `EXCLUDE_PATTERNS = */detail/*` (Phase 8's `diff/detail/` and any future internal headers stay out).
  - `EXTRACT_ALL = NO`, `EXTRACT_PRIVATE = NO`, `EXTRACT_STATIC = NO` — public Doxygen-block symbols only.
  - `GENERATE_LATEX = YES`, `LATEX_OUTPUT = latex`, `PDF_HYPERLINKS = YES`, `USE_PDFLATEX = YES`.
  - `GENERATE_HTML = YES` (for local browsing; not published).
  - `HAVE_DOT = YES`, `CLASS_DIAGRAMS = YES`, `COLLABORATION_GRAPH = YES`, `INCLUDE_GRAPH = YES`, `INCLUDED_BY_GRAPH = YES`, `CALL_GRAPH = NO`.
  - `WARN_IF_UNDOCUMENTED = YES`, `WARN_AS_ERROR = YES` (catches public symbols with no Doxygen block at all — paired with the `\since` lint in Group 6).
  - `OUTPUT_DIRECTORY = build/docs/doxygen` (relative to the source tree; matches the LaTeX `\input{}` path).
  - `PROJECT_NAME = "Grid Calculus"`, `PROJECT_NUMBER = @PROJECT_VERSION@` (CMake-substituted), `PROJECT_BRIEF = "Differential operators and functional evaluation on periodic crystalline grids"`.
- Update root `CMakeLists.txt` to substitute `@PROJECT_VERSION@` into a generated `Doxyfile.configured` under `${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile`. Only when `GRIDCALC_BUILD_DOCS` is `ON` (new build option, default `OFF` so contributors who haven't installed TeX Live or Doxygen aren't blocked).
- Add CMake `find_package(Doxygen)` and a `gridcalc-docs-doxygen` custom target that runs Doxygen on the configured Doxyfile.

**Commit message:** `build(docs): add Doxyfile and CMake docs target`

## Group 2 — LaTeX skeletons for both books

- Create `docs/common/preamble.tex` — shared package definitions (`amsmath`, `amssymb`, `siunitx`, `listings`, `xcolor`, `hyperref` with conservative settings, `tcolorbox` for callouts, project-local math macros).
- Create `docs/user-guide/main.tex` — `\documentclass{memoir}` skeleton:
  - Title page, copyright page, TOC.
  - `\input{../common/preamble.tex}`.
  - 12 chapter `\input{}`s under `docs/user-guide/chapters/`:
    - `01-hello-grid.tex` — Hello grid tutorial (Group 5).
    - `02-laplacian.tex` through `10-spectral-verification.tex` — Phase 1 through Phase 9 absorbed user-guide notes (Groups 3, 4).
    - `11-cahn-hilliard.tex`, `12-diamond-lattice.tex` — placeholder stubs for Phase 12 / 16.
  - Back matter: index (`\printindex`), references (placeholder).
- Create `docs/developer-note/main.tex` — `\documentclass{book}` skeleton:
  - Title page, copyright page, TOC.
  - `\input{../common/preamble.tex}`.
  - 10 chapter `\input{}`s under `docs/developer-note/chapters/`:
    - `00-overview.tex` — short orientation / how-to-read-this-book.
    - `01-laplacian.tex` through `09-spectral-verification.tex` — Phase 1 through Phase 9 absorbed developer-note notes (Group 4).
  - Back matter: `\input{build/docs/doxygen/latex/refman.tex}` for the auto-generated API reference (handled by Group 6 — Doxygen runs first in the docs build pipeline, then `latexmk` runs over the Developer Note).
- Create empty placeholder files for every chapter `\input{}` (`\chapter{Title}` + a TODO marker), so `latexmk` succeeds end-to-end before any content migration. Also create empty `docs/{user-guide,developer-note}/chapters/` directories if needed.
- The CMake `gridcalc-docs-userguide` and `gridcalc-docs-developernote` custom targets run `latexmk -pdf -output-directory=build/docs main.tex` from the respective `docs/<book>` directory. Build artifacts go under `build/docs/<book>/`.

**Commit message:** `feat(docs): add memoir + book LaTeX skeletons for both books`

## Group 3 — Migrate User Guide markdown notes to LaTeX chapters

- For each of the 9 markdown notes under `docs/user-guide/notes/phase-N-*.md`:
  1. Run `pandoc <note>.md -o docs/user-guide/chapters/<NN>-<slug>.tex --listings -f markdown -t latex` (the `--listings` flag uses the `listings` package for code blocks, matching the preamble).
  2. Hand-polish the generated `.tex`:
     - Replace pandoc's verbatim emoji / Unicode chars with LaTeX equivalents (`\psi`, `\nabla`, `\partial`, etc.).
     - Replace pandoc-generated `\hyperref` macros with project-local cross-references (`\Cref{chap:<slug>}`, `\Cref{sec:...}`).
     - Restructure code listings with consistent style (`lstset` defaults from `preamble.tex`).
     - Preserve original section structure (the existing markdown notes are already structured for absorption).
  3. Each chapter starts with a `\chapter{<Title>}` and uses `\section`/`\subsection` for the original `##` / `###` headings.
- Verify each chapter compiles cleanly under `latexmk -pdf` against the User Guide skeleton.
- Empty `docs/user-guide/notes/` directory once all 9 notes are migrated.
- Update `docs/user-guide/notes/README.md` (the explainer for the directory) to redirect future doc-content additions to the LaTeX chapters and remove the old worked-example skeleton.

**Commit message:** `docs(user-guide): absorb 9 phase notes into LaTeX chapters`

## Group 4 — Migrate Developer Note markdown notes to LaTeX chapters

- Same workflow as Group 3, applied to the 9 Developer Note markdown files under `docs/developer-note/notes/phase-N-*.md`.
- The Developer Note chapters are denser (five-section structure: Theory · Math derivation · Algorithm · Design decisions · References). Hand-polish is heavier here — math-mode display equations, citation re-formatting, file-path references with `\path{}` instead of pandoc's verbatim format.
- Empty `docs/developer-note/notes/` directory.
- Update `docs/developer-note/notes/README.md` similarly to Group 3.

**Commit message:** `docs(developer-note): absorb 9 phase notes into LaTeX chapters`

## Group 5 — "Hello grid" tutorial chapter + `docs/README.md`

- Author `docs/user-guide/chapters/01-hello-grid.tex` from scratch — the User Guide's lead chapter:
  - Section 1: "Setting up a grid and field" — `Grid`, `Field<double>`, sampler-lambda constructor (Phase 1).
  - Section 2: "First derivatives" — `diff::laplacian`, `diff::gradient`, `diff::divergence` at default `Order = 2` (Phases 1, 2).
  - Section 3: "Higher accuracy" — `<Order = 4>` (Phase 7); cost vs accuracy table.
  - Section 4: "Functionals" — `func::integrate`, `func::evaluate` with the Ginzburg–Landau worked example (Phases 3, 4).
  - Section 5: "Time integration" — `solve::diffuse` with the `RK4` tag (Phases 5, 6); a small heat-equation simulation.
  - Section 6: "Higher-order spatial operators" — `diff::biharmonic` and one mixed partial like `diff::d4dx2dy2<4>` (Phase 8).
  - Section 7: "Verification" — the FFT cross-check pattern using `spectral::laplacian` (Phase 9).
  - Each section is a few paragraphs + one short code listing. The chapter reads as a single continuous walkthrough; absorbed phase chapters (02–10) provide the per-feature reference depth.
- Create `docs/README.md` documenting:
  - Local build commands for each artifact (Doxygen, User Guide PDF, Developer Note PDF).
  - The `cmake --preset gcc-debug -DGRIDCALC_BUILD_DOCS=ON && cmake --build build/<preset> --target gridcalc-docs-all` flow.
  - The `scripts/build-docs.sh` convenience driver (created next bullet).
  - Required local toolchain: TeX Live (recommended + extra + fonts-recommended + pictures), latexmk, Doxygen ≥ 1.9, Graphviz.
  - Pointer to `scripts/get-docs.sh --ci` for fetching the latest CI-built PDFs without local toolchain.
- Create `scripts/build-docs.sh` — driver that does `cmake --build build/<preset> --target gridcalc-docs-doxygen gridcalc-docs-userguide gridcalc-docs-developernote` (or runs Doxygen + `latexmk` directly when called outside CMake). Replaces the deleted `scripts/render-docs.sh` for the local-render workflow.
- **Delete `scripts/render-docs.sh`** — retired this phase (per `requirements.md`).
- Update `scripts/get-docs.sh` to reflect the new artifact bundle (still named `gridcalc-docs-pdfs`; bundle now contains three artifacts — Doxygen output tarball, `user-guide.pdf`, `developer-note.pdf`).

**Commit message:** `docs: add Hello grid chapter, README, build-docs.sh; retire render-docs.sh`

## Group 6 — CI: full doc-build job + `\since` lint

- Update `.github/workflows/ci.yml`:
  - **Replace** the existing `Render Doc PDFs` job with a new `docs-build` job that:
    - Runs on `ubuntu-latest`.
    - Installs `texlive-latex-recommended texlive-latex-extra texlive-fonts-recommended texlive-pictures latexmk doxygen graphviz`.
    - Runs `cmake -B build/docs -DGRIDCALC_BUILD_DOCS=ON -DGRIDCALC_BUILD_TESTS=OFF`.
    - Runs `cmake --build build/docs --target gridcalc-docs-doxygen gridcalc-docs-userguide gridcalc-docs-developernote`.
    - Uploads three artifacts in one bundle named `gridcalc-docs-pdfs`: `build/docs/doxygen/` (tarball), `build/docs/user-guide/user-guide.pdf`, `build/docs/developer-note/developer-note.pdf`.
  - Add a new `docs-lint` job that:
    - Runs on `ubuntu-latest`, no toolchain install needed beyond Python 3.
    - Runs `python3 scripts/check-since.py include/gridcalc/`.
    - Also runs Doxygen with `WARN_AS_ERROR=YES` against the configured Doxyfile (catches missing-block symbols at compile time of the docs).
- Add `scripts/check-since.py`:
  - Walks `include/gridcalc/**/*.hpp` (excluding `detail/` subdirs).
  - For each file, parses public-API declarations preceded by a `///` Doxygen comment block. Identifies declarations matching: `template<...>`, `struct <Name>`, `class <Name>`, function declarations (return type + name + signature), `using <Name> = ...`, `static constexpr ... <Name>`.
  - Asserts the preceding Doxygen block contains an explicit `\since` line. Prints the file:line of any miss and exits with code 1 if any miss is found.
  - Includes a small allowlist for known un-tagged symbols (initially empty; filled if the lint surfaces any false positives during Phase 10's audit).
- Run both checks locally during Phase 10's audit step. Add any missing `\since` tags to the affected public headers before Group 6 commits.

**Commit message:** `ci(docs): replace render-docs job with full LaTeX build and add \since lint`

## Group 7 — Bookkeeping

- Bump `CMakeLists.txt` `project(gridcalc VERSION 0.9.0)` → `0.10.0`.
- Update `test/version_test.cpp` (`ReturnsZeroNineZero` → `ReturnsZeroTenZero`).
- `CHANGELOG.md`: new `## 0.10.0 — Phase 10: Documentation infrastructure` entry covering the Doxygen + LaTeX skeletons + chapter migration + retired pandoc render + new CI doc-build job + `\since` lint.
- `specs/STATUS.md`: bump "Last updated"; move Phase 10 row to Done; update "Where the project stands" + "Next action" (Phase 11 — Functional with up-to-4th-order gradients). Add new "Decisions worth knowing" entries: full LaTeX toolchain landed; `notes/` directories emptied; `\since` lint is a CI gate from Phase 10 forward.
- Update `specs/tech-stack.md` Documentation section to remove pandoc from the required toolchain (was added post-Phase-4 as the lightweight render's dep) and to note the `\since` lint policy.

**Commit message:** `chore: bump version to 0.10.0 and refresh STATUS for Phase 10`
