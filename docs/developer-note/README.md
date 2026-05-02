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

where `<slug>` matches the phase's spec-directory slug. The file follows a
**fixed five-section structure** (per `specs/CLAUDE.md` step 4d). Depth scales
with the phase's mathematical content; every section must be addressed (or
explicitly marked `N/A` with a one-line reason).

### Required sections

1. **Theory** — the mathematical objects and identities the code rests on
   (continuous PDE, function space, conservation law, lattice symmetry, …).
2. **Math derivation** — explicit derivation of the discrete formulas:
   Taylor-expansion / moment-matching for stencil weights, truncation-error
   analysis, eigenvalue spectra for trig manufactured solutions, stability
   arguments, resulting convergence order. Use display math (`$$ … $$`); the
   markdown will later be inlined into the LaTeX Developer Note where the
   math renders natively.
3. **Algorithm** — the procedure the code actually executes, with
   file-path references to the relevant headers (`include/gridcalc/...`) and
   tests (`test/...`). Note algorithmic complexity, data-layout decisions,
   non-trivial control-flow.
4. **Design decisions** — the non-obvious decisions taken in
   `requirements.md` distilled for posterity (one short paragraph each),
   including trade-offs considered and rejected. Cite `requirements.md`
   rather than copying it.
5. **References** — every theorem, formula, or algorithmic choice cited
   above gets an inline marker that resolves here. **At least one external
   reference per note.** Acceptable forms:
   - Journal / preprint paper: full bibliographic entry plus a permanent
     identifier (DOI URL `https://doi.org/...` preferred; arXiv ID
     acceptable for preprints).
   - Textbook: author, title, edition, publisher, year, plus the specific
     chapter / section / page being cited.
   - Permanent URL: a stable web reference (project documentation, official
     spec). No temporary links, social-media posts, or paywalled-without-DOI
     sources.

   Cross-references to other Developer Note phases or to the project's
   `specs/` files do not count toward the "at least one external" minimum.

### Skeleton

```markdown
# Developer Note — Phase N: <Phase title>

> Phase metadata. Spec dir: `specs/YYYY-MM-DD-phase-N-<slug>/`. PR: #<num>.
> Merged: YYYY-MM-DD.

## 1. Theory

<continuous setting; the mathematical object being approximated>

## 2. Math derivation

<explicit derivation; display math; truncation-error analysis; convergence
order>

$$
\partial_x f(x_i) \;\approx\; \frac{f(x_{i+1}) - f(x_{i-1})}{2 h}
\;+\; \mathcal{O}(h^2)
$$

(Truncation derivation here, citing reference [1].)

## 3. Algorithm

The discrete operator is implemented in
[`include/gridcalc/diff/Gradient.hpp`](../../../include/gridcalc/diff/Gradient.hpp).
Per grid point, one stencil application along each of the three axes; total
$O(N_x N_y N_z)$ work; periodic wrap is delegated to the input `Field`'s
`Policy`. Tested in
[`test/gradient_test.cpp`](../../../test/gradient_test.cpp).

## 4. Design decisions

- **Decision A.** One-paragraph summary citing
  [requirements.md](../../../specs/YYYY-MM-DD-phase-N-<slug>/requirements.md).

## 5. References

1. LeVeque, R. J. (2007). *Finite Difference Methods for Ordinary and Partial
   Differential Equations*. SIAM. <https://doi.org/10.1137/1.9780898717839>,
   §1.2 (centered differences).
2. <https://en.wikipedia.org/wiki/Finite_difference_coefficient> (permanent
   URL — coefficient tables for higher-order stencils).
```

## What does NOT go here

- User-facing API summaries — those go in `docs/user-guide/notes/`.
- Throwaway debugging notes or stale plans — the Developer Note is
  permanent; if a decision changes later, edit the existing note rather
  than appending.
- Citation-free assertions about correctness or performance. Every
  non-obvious claim needs a reference.

## Lifecycle

When Phase 10 lands, every existing `notes/phase-N-*.md` is absorbed into the
corresponding chapter of `docs/developer-note/<developer-note>.tex` and this
`notes/` directory is emptied. From Phase 11 forward, doc deliverables go
straight to the LaTeX sources; this placeholder tree is then retired.

See `specs/CLAUDE.md` step 4d (the doc-notes description) and step 4c (the
validation checklist) for the full process.
