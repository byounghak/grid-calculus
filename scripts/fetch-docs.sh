#!/usr/bin/env bash
#
# scripts/fetch-docs.sh
#
# Downloads the latest gridcalc-docs-pdfs artifact (user-guide.pdf,
# developer-note.pdf) from a successful CI run on the
# byounghak/grid-calculus repository and extracts it into build/docs/
# by default. Pairs with scripts/render-docs.sh, which produces the
# same files locally if pandoc + xelatex are installed.
#
# Usage:
#   ./scripts/fetch-docs.sh                       # latest successful main run
#   ./scripts/fetch-docs.sh --pr <num>            # latest run on the named PR
#   ./scripts/fetch-docs.sh --run <id>            # specific workflow run ID
#   ./scripts/fetch-docs.sh --out <dir>           # extract somewhere else
#   ./scripts/fetch-docs.sh --pr 7 --out ~/Desktop
#
# Requirements:
#   - gh (GitHub CLI), authenticated to the byounghak/grid-calculus repo.
#     Install: `brew install gh && gh auth login`.

set -euo pipefail

ARTIFACT_NAME="gridcalc-docs-pdfs"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="$ROOT/build/docs"
RUN_ID=""
PR_NUM=""

usage() {
  cat <<EOF
Usage: $(basename "$0") [--pr <num>] [--run <id>] [--out <dir>]

With no flags, downloads the latest successful CI run on the 'main'
branch and extracts the gridcalc-docs-pdfs artifact into build/docs/.

  --pr <num>     Use the latest successful run on the named PR's branch.
  --run <id>     Use a specific workflow run ID (takes priority over --pr).
  --out <dir>    Extract the artifacts here. Default: build/docs/.
  -h, --help     Show this help and exit.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pr)  PR_NUM="${2:-}"; shift 2 ;;
    --run) RUN_ID="${2:-}"; shift 2 ;;
    --out) OUT_DIR="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if ! command -v gh >/dev/null 2>&1; then
  echo "error: gh (GitHub CLI) not found on PATH." >&2
  echo "       install with: brew install gh && gh auth login" >&2
  exit 1
fi

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

run_id="$(resolve_run_id || true)"
if [[ -z "$run_id" ]]; then
  echo "error: no successful CI run found for the requested target." >&2
  exit 1
fi

mkdir -p "$OUT_DIR"
echo "fetching artifact '$ARTIFACT_NAME' from run $run_id"
echo "into $OUT_DIR ..."
gh run download "$run_id" --name "$ARTIFACT_NAME" --dir "$OUT_DIR"

echo "done. PDFs in $OUT_DIR:"
ls -1 "$OUT_DIR"/*.pdf 2>/dev/null || echo "  (no .pdf files extracted; check the artifact contents)"
