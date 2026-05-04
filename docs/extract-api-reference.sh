#!/usr/bin/env sh
# docs/extract-api-reference.sh
#
# Strips the standalone-document framing from Doxygen's refman.tex so its
# content portion can be \input{}'d into the Developer Note's appendix.
#
# Refman structure (Doxygen 1.16):
#   ...preamble...
#   \begin{document}
#   ...title page, \tableofcontents...
#   \chapter{Hierarchical Index}        <-- start of content
#   ...chapter / \input{} list...
#   \printindex
#   \end{document}
#
# This script extracts everything from the first `\chapter{Hierarchical
# Index}` line through the last `\input{}` directive before `\printindex`.

set -eu

SRC="${1:?usage: $0 <refman.tex> <api-reference.tex>}"
DEST="${2:?usage: $0 <refman.tex> <api-reference.tex>}"

if [ ! -f "$SRC" ]; then
    echo "extract-api-reference.sh: source file not found: $SRC" >&2
    exit 1
fi

# Use sed addresses: print from the line matching "\chapter{Hierarchical Index}"
# up to (but not including) the line matching "\printindex". The Doxygen
# preamble and \begin{document} / \tableofcontents are dropped. The trailing
# \end{document} is also dropped.
sed -n '/\\chapter{Hierarchical Index}/,/\\printindex/{ /\\printindex/!p; }' "$SRC" > "$DEST"

# Ensure the file is non-empty and looks valid.
if ! [ -s "$DEST" ]; then
    echo "extract-api-reference.sh: extracted file is empty: $DEST" >&2
    exit 1
fi
