# User Guide — placeholder markdown notes

This directory and its `notes/` subdirectory hold the **placeholder form** of
the User Guide until Phase 10 stands up the LaTeX skeleton (`memoir` class,
built with `latexmk -pdf` per `specs/tech-stack.md`).

## What goes here

From Phase 3 onward, every feature PR adds a single markdown file at:

```
docs/user-guide/notes/phase-N-<slug>.md
```

where `<slug>` matches the phase's spec-directory slug (e.g.,
`phase-2-gradient-divergence`). Each note targets the **production user**
driving the library and covers:

- What new public API was added or changed.
- Function signatures with brief docstrings (link to the headers; do not
  paraphrase the Doxygen blocks).
- A short, runnable worked example that exercises the new feature against a
  trig manufactured solution or similar.
- Any user-facing constraints (boundary policy, periodicity assumptions, grid
  layout, units).

Keep notes tight. The bar is "a user could read this and call the new API
correctly," not "a user could re-derive the math."

## What does NOT go here

- Internal design notes, math derivations, decision context — those go in
  `docs/developer-note/notes/`.
- Doxygen-style API tables — Phase 10 onward, the User Guide does not
  duplicate the API reference.

## Lifecycle

When Phase 10 lands, every existing `notes/phase-N-*.md` is absorbed into the
corresponding chapter of `docs/user-guide/<user-guide>.tex` and this `notes/`
directory is emptied. From Phase 11 forward, doc deliverables go straight to
the LaTeX sources; this placeholder tree is then retired.

See `specs/CLAUDE.md` step 4 (the doc-notes paragraph above 4b) and step 4c
(the validation checklist) for the full process.
