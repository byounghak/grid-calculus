# Developer Note — placeholder markdown notes

This directory and its `notes/` subdirectory hold the **placeholder form** of
the Developer Note until Phase 10 stands up the LaTeX skeleton (`book` class,
`\input{}`s Doxygen-generated subfiles, built with `latexmk -pdf` per
`specs/tech-stack.md`).

## What goes here

From Phase 3 onward, every feature PR adds a single markdown file at:

```
docs/developer-note/notes/phase-N-<slug>.md
```

where `<slug>` matches the phase's spec-directory slug. Each note targets an
**engineer reading the codebase later** and covers:

- The non-obvious design decisions taken in `requirements.md` for this phase
  (one short paragraph per decision; refer to the spec rather than copying).
- Math derivations and correctness arguments that don't fit in the Doxygen
  comments — e.g., truncation-error analyses, eigenvalue calculations,
  proofs of convergence order.
- Trade-offs considered and rejected, with enough context that a future
  contributor can re-evaluate them.
- Internal-only notes: invariants, data-layout decisions, performance
  hypotheses to verify in Phase 20.

Cite the relevant headers (`include/gridcalc/...`) and tests (`test/...`) by
path; do not duplicate the code or test bodies in the note.

## What does NOT go here

- User-facing API summaries — those go in `docs/user-guide/notes/`.
- Throwaway debugging notes or stale plans — the Developer Note is permanent;
  if a decision changes later, edit the existing note rather than appending.

## Lifecycle

When Phase 10 lands, every existing `notes/phase-N-*.md` is absorbed into the
corresponding chapter of `docs/developer-note/<developer-note>.tex` and this
`notes/` directory is emptied. From Phase 11 forward, doc deliverables go
straight to the LaTeX sources; this placeholder tree is then retired.

See `specs/CLAUDE.md` step 4 (the doc-notes paragraph above 4b) and step 4c
(the validation checklist) for the full process.
