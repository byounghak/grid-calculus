#!/usr/bin/env python3
"""\since-tag lint for public gridcalc headers.

Fails if any public-API symbol has a Doxygen `///` comment block but no
`\since` tag inside the block. Run from the repo root as:

    python3 scripts/check-since.py include/gridcalc/

The script walks every `.hpp` file under the given roots, skipping any
path that contains a `/detail/` segment (internal helpers; Phase 8's
`diff/detail/MultiIndexDerivative.hpp` is the canonical example).

For each header it looks at every Doxygen comment block and asks: does
the next non-blank, non-comment line below it look like a *declaration*
that a public consumer would care about?  If so, the block must contain
the substring `\\since` (or `@since`).  Each violation prints the
file:line of the declaration and contributes to a non-zero exit code.

The script is deliberately *narrow* on what it considers a public
declaration so that benign things like internal helpers, namespace
braces, and partial-template `;`-only forward declarations don't fire.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


# Lines starting with these patterns count as "looks like a public-API
# declaration" — i.e., a class/struct, free function (return type +
# name + signature), template, type alias, or namespace-scoped variable.
_PUBLIC_DECL = re.compile(
    r"""^\s*(
        template\s*<                # template <...> declaration
      | (?:inline\s+|constexpr\s+|static\s+|extern\s+)*
        (?:struct|class|enum)\s+\w+ # struct/class/enum
      | (?:inline\s+|constexpr\s+|static\s+|extern\s+)*
        (?:[\w:<>,\s\&\*]+?\s+)+    # return type (very permissive)
        \w+\s*\(                    # function name + open paren
      | using\s+\w+\s*=             # type alias
    )
    """,
    re.VERBOSE,
)

# Doxygen comment block opener.  We require triple-slash `///` (the
# project's house style); `/** */` and `//!` blocks are not used in
# gridcalc and are not scanned.
_DOXY_LINE = re.compile(r"^\s*///")


def is_under_detail(path: Path) -> bool:
    return any(part == "detail" for part in path.parts)


def check_file(path: Path) -> list[tuple[int, str]]:
    """Return list of (line_no, snippet) where a Doxygen block lacks \since."""
    violations: list[tuple[int, str]] = []
    lines = path.read_text(encoding="utf-8").splitlines()
    n = len(lines)
    i = 0
    while i < n:
        if not _DOXY_LINE.match(lines[i]):
            i += 1
            continue
        # Found the first line of a Doxygen block; collect contiguous
        # /// lines.
        block_start = i
        block: list[str] = []
        while i < n and _DOXY_LINE.match(lines[i]):
            block.append(lines[i])
            i += 1
        # Find the next non-blank, non-comment line — the candidate
        # declaration.  Skip blank lines and `// ...` non-Doxygen comments.
        decl_idx = i
        while decl_idx < n:
            line = lines[decl_idx]
            stripped = line.strip()
            if not stripped:
                decl_idx += 1
                continue
            # Plain `//` (non-Doxygen) comments — skip.
            if stripped.startswith("//") and not stripped.startswith("///"):
                decl_idx += 1
                continue
            break
        if decl_idx >= n:
            continue
        candidate = lines[decl_idx]
        if not _PUBLIC_DECL.match(candidate):
            continue
        # Skip declarations clearly inside a private/protected section
        # (`private:` or `protected:` access specifier — gridcalc public
        # headers don't use these much, but let's be safe).  We simply
        # don't check; the project convention is to only Doxygen-block
        # public symbols, so a private symbol with a /// block is rare
        # but acceptable.
        block_text = "\n".join(block)
        if "\\since" in block_text or "@since" in block_text:
            continue
        violations.append(
            (block_start + 1, candidate.strip()[:80])
        )
    return violations


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "roots",
        nargs="+",
        type=Path,
        help="Directories to recurse into (e.g. include/gridcalc/).",
    )
    args = p.parse_args(argv)

    total_violations = 0
    for root in args.roots:
        for hpp in sorted(root.rglob("*.hpp")):
            if is_under_detail(hpp):
                continue
            for line_no, snippet in check_file(hpp):
                print(f"{hpp}:{line_no}: missing \\since tag — {snippet}")
                total_violations += 1

    if total_violations:
        print(
            f"\ncheck-since: {total_violations} public symbol(s) missing "
            f"\\since tag.",
            file=sys.stderr,
        )
        return 1

    print("check-since: all public symbols carry \\since tags.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
