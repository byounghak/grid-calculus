# Requirements: Phase 10 — Documentation infrastructure

## Roadmap reference

> ## Phase 10 — Documentation infrastructure
>
> - **Goal.** Stand up the full documentation toolchain so every later phase can land its own incremental docs (API reference, User Guide chapters, Developer Note sections) into a working build. This phase is plumbing — each subsequent feature phase remains responsible for the doc deliverables it already lists (Cahn–Hilliard tutorial in Phase 12, diamond-lattice example in Phase 16, etc.); they fold in here rather than waiting for the v1.0 polish at Phase 22.
> - **Replaces the lightweight pandoc-based render.** A `pandoc + pdflatex` render (`scripts/render-docs.sh` + the CI `render-docs` job, introduced post-Phase-4) currently produces interim aggregate PDFs from the markdown notes. Phase 10 retires that script and CI job in favor of the full LaTeX skeleton below.
> - **Deliverables.**
>   - `docs/Doxyfile` configured for the public headers; Graphviz `dot` enabled for class / include / call graphs.
>   - LaTeX skeletons for the **Developer Note** (`book` class; `\input{}`s Doxygen-generated subfiles for the API reference and Graphviz diagrams) and the **User Guide** (`memoir` class; fully hand-authored). Both built with `latexmk -pdf`.
>   - **Absorb the markdown placeholder notes** under `docs/user-guide/notes/` and `docs/developer-note/notes/` (accumulated from Phase 3 onward, one pair per phase) into the corresponding `.tex` chapters. The two `notes/` directories are emptied at the end of this phase; from Phase 11 onward, doc deliverables go straight to the LaTeX sources.
>   - A "Hello grid" tutorial chapter in the User Guide covering the operators delivered through Phase 9 (scalar periodic FD, integration, callable functional, explicit/RK4 diffusion, higher-order accuracy, biharmonic + mixed partials, FFT verification). This chapter assembles content largely from the absorbed markdown notes.
>   - `docs/README.md` documenting the local build commands for each artifact.
>   - CI job: installs `texlive-latex-recommended texlive-latex-extra texlive-fonts-recommended texlive-pictures latexmk doxygen graphviz`; produces three artifacts on every PR — Doxygen output, Developer Note PDF, User Guide PDF. PDFs are not committed.
>   - **`\since` policy enforcement.** A CI lint scans the public headers and fails the build if any public class, function, member variable, or type alias is missing a `\since` tag. All symbols delivered through Phase 9 carry tags before this phase merges.
> - **Acceptance.**
>   - All three documentation artifacts build cleanly in CI on every PR.
>   - The `\since` lint blocks PRs that introduce public symbols without the tag.
>   - Build-cost budget: full docs build adds < 90 s to a clean CI run.
>   - The User Guide PDF contains the "Hello grid" chapter; the Developer Note PDF contains an auto-built API reference for everything in Phases 0–9 with no broken cross-references.

## Decisions made for this feature

- **Note migration: pandoc first-pass + hand polish.** Run `pandoc <note>.md -o <chapter>.tex` per markdown note as a starting point, then hand-edit each generated chapter for typography, cross-references, table layout, and consistency with the surrounding LaTeX skeleton. Cuts the migration time roughly in half vs hand-rewriting; preserves the structure of the markdown notes (which is already organized for the absorbed-into-LaTeX use case). The pandoc-generated `.tex` is committed only after hand-polish — no auto-generated artifacts in the repo. The 9 markdown pairs (`phase-1-laplacian` through `phase-9-spectral-verification`) become 9 chapters per book.
- **Doxygen → LaTeX integration: input `refman.tex`.** Doxygen runs with `GENERATE_LATEX=YES` and emits a complete `latex/refman.tex` plus per-symbol `.tex` snippets under `build/docs/doxygen/latex/`. The Developer Note's back matter `\input{}`s `build/docs/doxygen/latex/refman.tex` to embed the full auto-generated API reference. This matches the roadmap text verbatim and is the standard Doxygen-LaTeX integration pattern. Per-symbol partial inputs were considered and rejected: more granular but requires a per-chapter input map that grows as the API grows.
- **`scripts/render-docs.sh` and the CI `render-docs` job retire on Phase 10 merge.** Both are deleted in this PR. `scripts/get-docs.sh` continues to work — it now fetches the new full-LaTeX artifact bundle (`gridcalc-docs-pdfs` keeps the same name; bundle now contains three artifacts: Doxygen output tarball, User Guide PDF, Developer Note PDF). Local-render workflow gets a new `scripts/build-docs.sh` documented in `docs/README.md` that drives `latexmk` + Doxygen for both books. Pandoc is dropped from the project's required toolchain.
- **`\since` lint: Doxygen `WARN_IF_UNDOCUMENTED` + `scripts/check-since.py`.** Two complementary checks, both required to pass in CI:
  1. Doxygen's `WARN_IF_UNDOCUMENTED=YES` and `WARN_AS_ERROR=YES` (or equivalent post-process check) catch public-API symbols with no Doxygen block at all.
  2. A small Python script `scripts/check-since.py` parses every public header (`include/gridcalc/**/*.hpp` excluding `detail/`), identifies declarations that should carry Doxygen blocks, and verifies each block contains an explicit `\since` line. The script is regex-based, scoped narrowly: every `template <...>`, `struct`, `class`, `void/auto/inline/... <name>(...)`, `using <name> =`, etc. preceded by a `///` Doxygen block must have `\since` somewhere in that block.
  3. Both checks run as a CI job named `docs-lint`. The job is fast (~5–10 s) and runs alongside the GCC / Clang / docs-build jobs.
- **LaTeX project layout.**
  - **User Guide** (`memoir` class, hand-authored): `docs/user-guide/main.tex` is the entry point; chapters live as `docs/user-guide/chapters/01-hello-grid.tex` through `12-spectral-verification.tex` (chapters 02–10 mirror Phase 1 through Phase 9 absorbed notes; chapter 01 is the "Hello grid" quickstart; chapters 11–12 are appendices for now). Front matter: title page, copyright, TOC. Back matter: index, license, references. Built with `latexmk -pdf`.
  - **Developer Note** (`book` class with `doxygen.sty` for the API-reference back matter): `docs/developer-note/main.tex` is the entry point; chapters live as `docs/developer-note/chapters/01-overview.tex` through `10-spectral-verification.tex` (chapters 01–10 mirror Phase 1–9 absorbed notes plus the new overview chapter); back matter `\input{}`s the Doxygen-generated `refman.tex`. Built with `latexmk -pdf`.
  - Both books share `docs/common/preamble.tex` for repeated package / macro definitions (math, code, links).
  - **`pdflatex` is locked** as the engine for both books (per `tech-stack.md`), even though the lightweight render used `xelatex`. The hand-authored LaTeX uses standard ASCII-mode math + named LaTeX commands (`\psi`, `\nabla`) for Greek and operators; no Unicode characters in the source.
- **"Hello grid" chapter scope.** Chapter 01 of the User Guide. Walks the reader through a single end-to-end example covering, in sequence: setting up a `Grid` and `Field<double>` (Phase 1); applying `laplacian` and `gradient` at default order 2 (Phases 1, 2); switching to `Order = 4` (Phase 7); evaluating a Ginzburg–Landau functional with `func::evaluate` (Phase 4); time-integrating diffusion via `solve::diffuse` with `RK4` (Phases 5, 6); a higher-order partial with `diff::d4dx2dy2<4>` and the biharmonic (Phase 8); and verifying against `spectral::laplacian` (Phase 9). The example is one continuous code listing built up incrementally, not a phase-by-phase summary. The phase-by-phase reference content lives in chapters 02–10.
- **Doxygen scope.** Doxygen is run over `include/gridcalc/` (public headers only); `internal` and `detail` namespaces are excluded via `EXCLUDE_PATTERNS = */detail/*`. Test code, examples, and benchmarks are excluded. Graphviz `dot` is enabled for `CLASS_DIAGRAMS=YES`, `COLLABORATION_GRAPH=YES`, `INCLUDE_GRAPH=YES`, `INCLUDED_BY_GRAPH=YES`, `CALL_GRAPH=NO` (call graphs add ~30s and aren't useful for a header-only library).
- **`\since` tag audit.** Before merging, every public symbol shipped through Phase 9 must carry a `\since` tag. Existing per-phase headers should already have them — the lint will catch any miss. The audit is part of the validation checklist.
- **Build-cost budget: < 90 s added to a clean CI run.** Phase 9's CI runtime is ~6 min on Ubuntu Clang (heavy by the FD–FFT cross-check fixture). Phase 10's docs build (Doxygen LaTeX + 2 × `latexmk` runs) is expected to add ~30–60 s wall-clock for the docs job. The `\since` lint adds ~5–10 s. Total well under the < 90 s budget. The new docs job runs in parallel with the existing GCC / Clang test jobs in CI.
- **Version bump: `0.10.0`.** No public-API change; new build infrastructure + tooling. SemVer minor bump.

## Non-goals

- **No content rewrite.** Existing markdown notes are migrated as-is structurally; hand-polish is for typography and cross-references, not factual content changes.
- **No new feature documentation.** Feature work belongs in feature phases; Phase 10 is plumbing only.
- **No HTML docs site.** The repository is private; PDFs are the only narrative-doc output. Doxygen produces HTML for local browsing if desired but is not published.
- **No Sphinx.** Decided in Phase 4; not revisited here.
- **No Doxygen call-graphs.** Adds ~30 s to the docs build for little reader benefit on a header-only library.
- **No XeLaTeX migration.** `pdflatex` is locked per `tech-stack.md`. The lightweight render's XeLaTeX (chosen for Unicode-tolerant input) is replaced by hand-authored LaTeX that uses LaTeX commands for non-ASCII glyphs.
- **No biblatex / hyperref deep customization.** Standard LaTeX cross-references; the absorbed notes' inline citations carry over with minor formatting edits.

## Deferred questions

- **Public docs hosting.** Repo is private; not publishing. Decision deferred indefinitely.
- **Cahn–Hilliard tutorial chapter.** Phase 12 owns the CH demo and its accompanying tutorial chapter; Phase 10 leaves chapter 11 of the User Guide as a stub placeholder for it.
- **Diamond-lattice example chapter.** Phase 16. Same pattern.
- **`\since` lint on private-namespace declarations.** Today's lint scans public headers only. Internal `detail::` declarations stay un-tagged for now.
- **PDF release-asset publication.** Phase 22 — the v1.0 release. Until then, PDFs are CI build artifacts only.
