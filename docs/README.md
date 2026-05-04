# Documentation builds

This directory holds the Phase 10 LaTeX skeletons for the **User Guide**
(`memoir`) and the **Developer Note** (`book`), the Doxygen
configuration, and the per-chapter content (most of which was migrated
from the per-phase markdown notes that lived under
`docs/{user-guide,developer-note}/notes/` through Phase 9).

## What gets built

Three artifacts, all under `build/<preset>/docs/` after a docs build:

- `doxygen/html/` — Doxygen HTML output (browse locally; not published).
- `doxygen/latex/` — Doxygen LaTeX output (`refman.tex` and per-symbol
  `.tex` files). The Developer Note `\input{}`s the chapter content
  from this tree at build time.
- `user-guide/main.pdf` — User Guide PDF.
- `developer-note/main.pdf` — Developer Note PDF (with the auto-built
  API reference appendix).

## Required local toolchain

- A C++17 compiler (CMake configure-time only — docs don't compile any
  C++ source).
- CMake 3.20+.
- TeX Live with `texlive-latex-recommended`, `texlive-latex-extra`,
  `texlive-fonts-recommended`, `texlive-pictures` packages installed.
  On macOS: `brew install --cask mactex-no-gui`.
- `latexmk` (typically bundled with TeX Live).
- Doxygen ≥ 1.9. On macOS: `brew install doxygen`.
- Graphviz for class / include diagrams. On macOS: `brew install graphviz`.

## Local build

From the repo root:

```sh
# Option A: one-shot driver (no CMake invocation needed).
./scripts/build-docs.sh

# Option B: via CMake explicitly.
cmake -B build/docs -DGRIDCALC_BUILD_DOCS=ON -DGRIDCALC_BUILD_TESTS=OFF
cmake --build build/docs --target gridcalc-docs-all
```

Both options produce all three artifacts under `build/docs/docs/`.

## CI build

GitHub Actions runs the same targets on every PR. The
`gridcalc-docs-pdfs` workflow artifact bundles the three outputs. To
download the latest artifact without a local toolchain:

```sh
./scripts/get-docs.sh --ci          # latest run on main
./scripts/get-docs.sh --pr 12       # latest run on PR #12
```

## When to add content

- **New public symbol** → put a Doxygen block on the declaration. The
  `docs-lint` CI job (Phase 10 onward) fails the build if any public
  symbol is missing a `\since` tag.
- **New User Guide chapter** → drop a new `docs/user-guide/chapters/NN-slug.tex`
  and add `\input{}` it from `docs/user-guide/main.tex`.
- **New Developer Note chapter** → same pattern under
  `docs/developer-note/chapters/NN-slug.tex`.
- The pre-Phase-10 placeholder `docs/{user-guide,developer-note}/notes/`
  directories are intentionally empty and exist only as guard rails for
  pull requests that try to add markdown notes — feature phases should
  add LaTeX chapters directly.

## Engine

`pdflatex` (per `specs/tech-stack.md`). The migrated chapters use
`inputenc utf8` + a `\DeclareUnicodeCharacter` allowlist in
`docs/common/preamble.tex` for the small set of Unicode glyphs the
markdown sources contained outside math mode. The Developer Note also
loads `hanging` with the `notquote` option to keep Doxygen's transitive
include from clobbering `f'(x)` in math expressions.
