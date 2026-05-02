# Claude.md — Feature workflow

When the user asks you to start the next phase of work (or otherwise begin a new feature in this repository), follow this workflow exactly. Every step happens before you start writing implementation code.

## 1. Find the next phase on the roadmap

Read [roadmap.md](roadmap.md) and identify the lowest-numbered phase whose work has not yet been delivered. Use the existing branch list, the git history, and the `specs/` directory contents to decide which phases are already done.

State the phase number and title to the user, and confirm before proceeding ("Next up is Phase 11 — Functional with up-to-4th-order gradients. Confirm?"). If the user wants to skip ahead or pick a different phase, follow their direction.

## 2. Open the branch first

Open a git branch named `YYYY-MM-DD-feature-name`, where:

- `YYYY-MM-DD` is today's date.
- `feature-name` is a short kebab-case slug derived from the phase title (examples: `phase-1-laplacian`, `phase-9-spectral-verification`, `phase-16-multi-atom-basis`).

The branch must exist **before** the clarifying questions in step 3, so that any commits you make during planning (e.g., adding the spec files themselves) land on the right branch.

## 3. Ask the user about the feature spec

Use the `AskUserQuestion` tool to gather the decisions needed to write `requirements.md`. Group the questions sensibly (one per major decision). Typical topics:

- Scope adjustments — anything in the roadmap-phase deliverables to drop, defer, or add?
- API choices not pinned by the roadmap or by [tech-stack.md](tech-stack.md) (signatures, naming, default arguments).
- Test data and reference inputs (which manufactured solution? which grid sizes for the convergence sweep?).
- Explicit out-of-scope items the user wants excluded.
- Anything the user has changed their mind on since the roadmap was written.

Do not skip this step even if the phase looks self-explanatory — there is almost always at least one open question worth pinning down.

## 4. Create the spec directory

Create `specs/YYYY-MM-DD-feature-name/` (same name as the branch). Inside it, write three files described below. All three may be revised during feature work; the bar is that they are accurate at merge time.

### 4a. `plan.md` — numbered, commit-sized task groups

Break the phase into a numbered list of groups, where each group is a single commit's worth of related work. A typical phase has 3–8 groups; a single PR contains all of them. Order groups so the build stays green between commits where practical (e.g., add a stub header before the first test that uses it).

Format:

```markdown
# Plan: <phase title>

## Group 1 — <short verb-phrase title>

- Task A
- Task B

**Commit message:** `<type>: <short imperative summary>`

## Group 2 — <short verb-phrase title>

- Task C
- ...
```

`<type>` follows conventional-commits-style prefixes appropriate to the change (`feat`, `test`, `docs`, `build`, `refactor`, `perf`, `chore`).

### 4b. `requirements.md` — scope, decisions, context (lightweight)

Do not write a full design doc. Capture only what is non-obvious from the roadmap or the global specs. Structure:

```markdown
# Requirements: <phase title>

## Roadmap reference

> Verbatim copy of the matching phase's Goal / Deliverables / Acceptance bullets from specs/roadmap.md, as a blockquote.

## Decisions made for this feature

- <each decision the user took during step 3, one bullet each>

## Non-goals

- <explicit out-of-scope items so future readers do not expand scope>

## Deferred questions

- <open questions punted to a later phase, named with the phase number that will address them>
```

Refer to [mission.md](mission.md) (audience, quality bar, in/out of scope) and [tech-stack.md](tech-stack.md) (language, dependencies, code style) for global context — cite them, do not restate.

### 4c. `validation.md` — merge-readiness checklist

A markdown tickbox checklist. Each item is a concrete, verifiable check that can be ticked during PR review. Group by category. Pull every Acceptance bullet from the roadmap phase into this list as its own check; add the standard build/test/docs items below.

Format:

```markdown
# Validation: <phase title>

## Build

- [ ] CI green on all blocking compilers (GCC, Clang, Apple Clang, MSVC).
- [ ] Warnings clean under `-Wall -Wextra -Wpedantic -Wconversion` (and MSVC equivalents).
- [ ] `clang-format` and `clang-tidy` pass.

## Tests

- [ ] <one bullet per Acceptance criterion from the roadmap phase>
- [ ] Convergence-order test passes within the documented tolerance (if applicable).

## Documentation

- [ ] Every new public symbol carries a Doxygen block including `\brief`, `\param`, `\returns`, `\since`.
- [ ] User Guide note added at `docs/user-guide/notes/phase-N-<slug>.md` covering user-facing API additions, with a worked example. Required from Phase 3 onward (placeholder markdown until Phase 10 absorbs it into the LaTeX User Guide).
- [ ] Developer Note added at `docs/developer-note/notes/phase-N-<slug>.md` covering the design decisions taken in `requirements.md` plus any internal-only notes, math derivations, or correctness arguments. Same placeholder convention as the User Guide note.
- [ ] `CHANGELOG.md` entry under the next minor version.

## Performance (if applicable)

- [ ] Benchmark numbers recorded under `bench/baselines/`.
- [ ] No regression beyond the documented threshold versus the prior baseline.
```

The PR does not merge with any unticked items.

### 4d. Doc-notes for User Guide and Developer Note (Phase 3 onward)

In addition to the three feature spec files in `specs/<phase>/`, every phase from Phase 3 onward adds a paired markdown note in the same PR as the implementation:

- `docs/user-guide/notes/phase-N-<slug>.md` — user-facing summary of the phase: what new public API was added, how to call it, what the inputs/outputs look like, with a short worked example. Audience is the production user driving the library.
- `docs/developer-note/notes/phase-N-<slug>.md` — implementation/design summary: the decisions taken in `requirements.md` distilled for an engineer reading the codebase later, plus any internal-only notes, math derivations, or correctness arguments.

The `<slug>` matches the phase's branch / spec-directory slug (e.g., `phase-2-gradient-divergence` → `phase-2-gradient-divergence.md`).

These markdown notes are the **placeholder** form. Phase 10 stands up the LaTeX skeletons for the User Guide (`memoir`) and the Developer Note (`book`), at which point the markdown notes are absorbed into the corresponding `.tex` chapters and the `notes/` directories are emptied. Until Phase 10 lands, markdown is the canonical narrative-doc format and the LaTeX toolchain is not yet built.

The two `notes/` directories carry their own `README.md` describing what does and does not belong in each. `validation.md` (4c) enforces the presence of both notes.

## 5. Reference, don't duplicate

The three project-wide spec files are the authoritative source of truth:

- [mission.md](mission.md) — purpose, audience, quality bar, in-scope / out-of-scope.
- [tech-stack.md](tech-stack.md) — language, dependencies, build, code style, documentation toolchain.
- [roadmap.md](roadmap.md) — phase sequence with Goal / Deliverables / Acceptance.

Feature specs cite these documents, they do not repeat them. If a feature requires a change to one of the global specs, edit that spec in the same branch and call out the change in `requirements.md` under "Decisions made for this feature."

## 6. After opening the PR — watch CI

Once the PR is open, watch the CI run on it (e.g. `gh pr checks <num> --watch`) until every blocking job reports a result. Do not declare the phase complete or recommend a merge until every blocking job for the current acceptance bar is green. If a job fails, investigate and fix on the same branch — do not bypass with `--no-verify`, `--admin` merges, or by relaxing the acceptance bar without an explicit `requirements.md` revision.
