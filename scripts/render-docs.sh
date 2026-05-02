#!/usr/bin/env bash
#
# scripts/render-docs.sh
#
# Renders the markdown User Guide and Developer Note notes into two PDFs:
#   build/docs/user-guide.pdf
#   build/docs/developer-note.pdf
#
# This is the **lightweight, pre-Phase-10 PDF render** introduced after
# Phase 4. Phase 10 replaces it with the full LaTeX skeleton (`memoir` for
# the User Guide, `book` for the Developer Note, driven by `latexmk`); see
# specs/roadmap.md Phase 10 and specs/tech-stack.md "Documentation".
#
# Requirements:
#   - pandoc                                            (markdown -> LaTeX -> PDF)
#   - xelatex (texlive-xetex on Ubuntu)                 (PDF engine; Unicode-native)
#
# Why xelatex rather than pdflatex: the doc-notes contain Greek letters and
# other Unicode glyphs both inside and outside math mode (e.g.,
# `O(ε log n)` in code spans). pdflatex requires explicit inputenc setup
# for non-ASCII; xelatex handles Unicode natively. Phase 10's full LaTeX
# skeleton stays on pdflatex (per specs/tech-stack.md), since it
# hand-authors all .tex sources and can quote Greek letters via LaTeX
# commands; this lightweight render is friendlier to whatever the author
# writes in markdown.
#
# CI installs both. Locally on macOS:
#   brew install pandoc
#   brew install --cask mactex-no-gui   # if xelatex is not already on PATH
#
# Output PDFs are gitignored (.gitignore already excludes *.pdf and build/).

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="$ROOT/build/docs"
mkdir -p "$OUT_DIR"

if ! command -v pandoc >/dev/null 2>&1; then
  echo "error: pandoc not found on PATH. See scripts/render-docs.sh header for install hints." >&2
  exit 1
fi
if ! command -v xelatex >/dev/null 2>&1; then
  echo "error: xelatex not found on PATH. See scripts/render-docs.sh header for install hints." >&2
  exit 1
fi

render_book() {
  local title="$1"
  local input_dir="$2"
  local output="$3"

  # Sort notes lexically; phase numbers are single-digit through Phase 9
  # (Phase 10 absorbs the markdown notes into LaTeX), so "phase-3-..." <
  # "phase-4-..." < ... naturally.
  local notes=()
  while IFS= read -r -d '' f; do
    notes+=("$f")
  done < <(find "$input_dir" -maxdepth 1 -name 'phase-*.md' -print0 | sort -z)

  if (( ${#notes[@]} == 0 )); then
    echo "no phase-*.md notes found in $input_dir; skipping $output"
    return
  fi

  echo "rendering $output from ${#notes[@]} note(s):"
  for f in "${notes[@]}"; do
    echo "  - $(basename "$f")"
  done

  pandoc \
    --metadata title="$title" \
    --metadata date="$(date -u +%Y-%m-%d)" \
    --toc \
    --toc-depth=2 \
    --number-sections \
    -V geometry:margin=1in \
    -V documentclass=report \
    -V colorlinks=true \
    -V linkcolor=blue \
    -V urlcolor=blue \
    --pdf-engine=xelatex \
    -o "$output" \
    "${notes[@]}"
}

render_book "Grid Calculus — User Guide" \
  "$ROOT/docs/user-guide/notes" \
  "$OUT_DIR/user-guide.pdf"

render_book "Grid Calculus — Developer Note" \
  "$ROOT/docs/developer-note/notes" \
  "$OUT_DIR/developer-note.pdf"

echo "done. PDFs are in $OUT_DIR"
