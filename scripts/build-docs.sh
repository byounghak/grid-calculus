#!/usr/bin/env bash
#
# scripts/build-docs.sh
#
# Builds all three documentation artifacts via the Phase 10 CMake docs
# pipeline:
#   - build/docs/docs/doxygen/{html,latex}/
#   - build/docs/docs/user-guide/main.pdf
#   - build/docs/docs/developer-note/main.pdf
#
# Replaces scripts/render-docs.sh (retired in Phase 10 alongside the
# pandoc-based interim render).
#
# Usage:
#   ./scripts/build-docs.sh                # default preset name `docs`
#   ./scripts/build-docs.sh <preset-name>  # override the build directory name
#
# Requirements:
#   - cmake, ninja
#   - TeX Live (latexmk + recommended packages); see docs/README.md
#   - Doxygen, Graphviz dot

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PRESET="${1:-docs}"
BUILD_DIR="$ROOT/build/$PRESET"

cd "$ROOT"

if ! command -v cmake >/dev/null 2>&1; then
    echo "error: cmake not found on PATH." >&2
    exit 1
fi
if ! command -v doxygen >/dev/null 2>&1; then
    echo "error: doxygen not found on PATH. See docs/README.md for install hints." >&2
    exit 1
fi
if ! command -v latexmk >/dev/null 2>&1; then
    echo "error: latexmk not found on PATH. See docs/README.md for install hints." >&2
    exit 1
fi

cmake -B "$BUILD_DIR" \
      -DGRIDCALC_BUILD_DOCS=ON \
      -DGRIDCALC_BUILD_TESTS=OFF \
      -G Ninja
cmake --build "$BUILD_DIR" --target gridcalc-docs-all

echo
echo "PDFs in $BUILD_DIR/docs/:"
ls -la "$BUILD_DIR"/docs/{user-guide,developer-note}/main.pdf 2>/dev/null || true
