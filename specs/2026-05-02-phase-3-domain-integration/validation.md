# Validation: Phase 3 — Domain integration (∫ over the grid)

## Build

- [ ] CI green on Ubuntu GCC and Ubuntu Clang (Phase 21 widens to Apple Clang + MSVC).
- [ ] Warnings clean under `-Wall -Wextra -Wpedantic -Wconversion`.
- [ ] `clang-format` and `clang-tidy` pass (subject to the codebase-wide hygiene gap noted in Phase 2's STATUS — no new regressions introduced).
- [ ] `gridcalc` target remains INTERFACE — no new source files.
- [ ] No new top-level dependency added (no OpenMP wiring at Phase 3).

## Tests

- [ ] `IntegrateTest.UnitFieldRecoversDomainVolume` passes — `∫ 1 dV` equals `Lx*Ly*Lz` to machine precision on an anisotropic, non-cubic Grid.
- [ ] `IntegrateTest.SinOverPeriodIsZero` passes — `∫ sin(kx) dV` over `[0, 2π]³` for integer `kx` is 0 within ≤ 1e-12.
- [ ] `IntegrateTest.ManufacturedSolutionPeriodicProduct` passes — `∫ cos(x) sin(2y) sin(3z) dV` over `[0, 2π]³` is 0 within machine precision.
- [ ] `IntegrateTest.PairwiseAndKahanAgree` passes — both reducers produce results agreeing within ≤ 1e-13 relative.
- [ ] `IntegrateTest.DeterministicAcrossInvocations` passes — back-to-back calls return bit-identical doubles. (Phase 20 replaces this with `OMP_NUM_THREADS`-varying invocations.)
- [ ] All Phase 1 + Phase 2 tests still pass.

## Documentation

- [ ] Every new public symbol (`func::integrate` overloads, `func::Pairwise`, `func::Kahan`) carries a Doxygen block with `\brief`, `\param`, `\returns`, `\since 0.3.0`.
- [ ] User Guide note added at `docs/user-guide/notes/phase-3-domain-integration.md` covering the signature, the two reducer tags, the periodic midpoint rule, and a worked example.
- [ ] Developer Note added at `docs/developer-note/notes/phase-3-domain-integration.md` with the five-section structure (Theory · Math derivation · Algorithm · Design decisions · References) per `specs/CLAUDE.md` step 4d. References section non-empty: at least one external citation (DOI URL, full bibliographic entry, or permanent URL).
- [ ] `CHANGELOG.md` entry under `0.3.0` lists: `func::integrate`, `Pairwise`, `Kahan`.
- [ ] `STATUS.md` updated: Phase 3 row marked Done; Phase 4 stated as the new Next action; "Last updated" bumped.

## Performance

- N/A at Phase 3. Performance work is centralized in Phase 20.
