#!/usr/bin/env bash
#
# scripts/get-docs.sh
#
# One-stop entry point for the gridcalc User Guide and Developer Note PDFs.
# Tries a local pandoc render first; falls back to downloading the
# gridcalc-docs-pdfs artifact from the most recent successful CI run on
# main if pandoc / xelatex are not installed.
#
# Usage:
#   ./scripts/get-docs.sh                     # auto: local if tools available, else CI fetch
#   ./scripts/get-docs.sh --local             # force local render (errors if tools missing)
#   ./scripts/get-docs.sh --ci                # skip local; fetch CI artifact
#   ./scripts/get-docs.sh --pr <num>          # implies --ci; from named PR's latest run
#   ./scripts/get-docs.sh --run <id>          # implies --ci; specific workflow run
#   ./scripts/get-docs.sh --out <dir>         # output directory (default: build/docs/)
#   ./scripts/get-docs.sh -h | --help         # show usage
#
# Requirements (one or the other; both is fine):
#   - Local render: pandoc + xelatex.
#       brew install pandoc
#       brew install --cask mactex-no-gui
#   - CI fetch: gh (GitHub CLI), authenticated to byounghak/grid-calculus.
#       brew install gh && gh auth login

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEFAULT_OUT_DIR="$ROOT/build/docs"
OUT_DIR="$DEFAULT_OUT_DIR"
ARTIFACT_NAME="gridcalc-docs-pdfs"
MODE="auto"
RUN_ID=""
PR_NUM=""

usage() {
  cat <<EOF
Usage: $(basename "$0") [--local|--ci] [--pr <num>] [--run <id>] [--out <dir>]

With no flags, tries a local pandoc render first; falls back to the latest
successful CI run on 'main' if pandoc/xelatex are not installed locally.

  --local        Force local render. Errors if pandoc or xelatex is missing.
  --ci           Skip the local-render attempt; fetch from CI directly.
  --pr <num>     Implies --ci. Fetch the named PR's latest successful run.
  --run <id>     Implies --ci. Fetch a specific workflow run by ID.
  --out <dir>    Output directory (default: build/docs/).
  -h, --help     Show this help and exit.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --local) MODE="local"; shift ;;
    --ci)    MODE="ci"; shift ;;
    --pr)    PR_NUM="${2:-}"; MODE="ci"; shift 2 ;;
    --run)   RUN_ID="${2:-}"; MODE="ci"; shift 2 ;;
    --out)   OUT_DIR="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

have_local_render() {
  command -v pandoc >/dev/null 2>&1 && command -v xelatex >/dev/null 2>&1
}

run_local_render() {
  echo "running local pandoc render"
  # render-docs.sh writes to build/docs/. Mirror to OUT_DIR if different.
  "$ROOT/scripts/render-docs.sh"
  if [[ "$OUT_DIR" != "$DEFAULT_OUT_DIR" ]]; then
    mkdir -p "$OUT_DIR"
    cp "$DEFAULT_OUT_DIR/user-guide.pdf" "$DEFAULT_OUT_DIR/developer-note.pdf" "$OUT_DIR/"
  fi
}

resolve_run_id() {
  if [[ -n "$RUN_ID" ]]; then
    echo "$RUN_ID"
    return
  fi
  if [[ -n "$PR_NUM" ]]; then
    local head_ref
    head_ref=$(gh pr view "$PR_NUM" --json headRefName --jq '.headRefName')
    if [[ -z "$head_ref" || "$head_ref" == "null" ]]; then
      echo "error: could not resolve head branch for PR #$PR_NUM" >&2
      exit 1
    fi
    gh run list --branch "$head_ref" --workflow ci.yml --status success \
      --limit 1 --json databaseId --jq '.[0].databaseId // empty'
    return
  fi
  gh run list --branch main --workflow ci.yml --status success \
    --limit 1 --json databaseId --jq '.[0].databaseId // empty'
}

run_ci_fetch() {
  if ! command -v gh >/dev/null 2>&1; then
    echo "error: gh (GitHub CLI) not found on PATH." >&2
    echo "       install with: brew install gh && gh auth login" >&2
    exit 1
  fi
  local run_id
  run_id="$(resolve_run_id || true)"
  if [[ -z "$run_id" ]]; then
    echo "error: no successful CI run found for the requested target." >&2
    exit 1
  fi
  mkdir -p "$OUT_DIR"
  echo "fetching artifact '$ARTIFACT_NAME' from run $run_id into $OUT_DIR ..."
  gh run download "$run_id" --name "$ARTIFACT_NAME" --dir "$OUT_DIR"
}

case "$MODE" in
  local)
    if ! have_local_render; then
      echo "error: --local requested but pandoc and/or xelatex not on PATH." >&2
      echo "       install with: brew install pandoc; brew install --cask mactex-no-gui" >&2
      exit 1
    fi
    run_local_render
    ;;
  ci)
    run_ci_fetch
    ;;
  auto)
    if have_local_render; then
      run_local_render
    else
      echo "pandoc/xelatex not found locally; falling back to CI artifact ..."
      run_ci_fetch
    fi
    ;;
esac

echo "done. PDFs in $OUT_DIR:"
ls -1 "$OUT_DIR"/*.pdf 2>/dev/null || echo "  (no .pdf files; check the steps above)"
