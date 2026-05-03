# Validation: Phase 9 — Spectral verification harness

## Build

- [ ] CI green on Ubuntu GCC and Ubuntu Clang with `GRIDCALC_USE_FFT=ON` (the default).
- [ ] CI also green with `GRIDCALC_USE_FFT=OFF` (skips spectral tests; FD tests pass).
- [ ] Warnings clean under `-Wall -Wextra -Wpedantic -Wconversion`. PocketFFT's vendored header is included with the `SYSTEM` modifier so its warnings do not fire on our code.
- [ ] `clang-format` and `clang-tidy` pass on the new project files. PocketFFT's vendored header is excluded from clang-tidy via `.clang-tidy` `HeaderFilterRegex`.
- [ ] `gridcalc` target remains INTERFACE — no new source files.
- [ ] `extern/pocketfft/VERSION.txt` records the upstream revision SHA + URL + snapshot date.

## Tests

- [ ] **Spectral basic / FFT round-trip.** `forwardR2C` followed by `backwardC2R` on $\sin(x)\sin(y)\sin(z)$ recovers the input within $5\times10^{-13}$ max-norm at $N = 32$.
- [ ] **Named wrappers smoke tests.**
  - `spectral::laplacian` on $\sin(\mathbf{k}\cdot\mathbf{x})$ for $\mathbf{k} = (1, 1, 1)$ recovers $-|\mathbf{k}|^2 \psi = -3 \psi$ within $5\times10^{-12}$ at $N = 32$.
  - `spectral::biharmonic` recovers $|\mathbf{k}|^4 \psi = 9 \psi$ within $5\times10^{-11}$ (looser than the laplacian bound because $|\mathbf{k}|^4 = 9$ amplifies round-off vs $|\mathbf{k}|^2 = 3$).
  - `spectral::gradient` recovers $\mathbf{k} \cos(\mathbf{k}\cdot\mathbf{x}) = (1, 1, 1) \cos(x+y+z)$ within $5\times10^{-12}$ component-wise.
- [ ] **FD–FFT cross-check fixture (the permanent CI gate).** For every (FD operator, Order) pair listed below, the log-log slope of the FD–FFT discrepancy on $\psi = \sin(x + y + z)$ over $N \in \{16, 32, 64, 128\}$ lies in `[Order − 0.5, Order + 0.5]`.
  - Phase 1 (2 cases): `diff::laplacian<{2, 4}>`.
  - Phase 2 (4 cases): `diff::gradient<{2, 4}>`, `diff::divergence<{2, 4}>`.
  - Phase 8 (64 cases): all 31 `d`-prefix partials × `Order ∈ {2, 4}` (62 cases) + `diff::biharmonic<{2, 4}>` (2 cases).
  - **Total: 70 cross-check tests.**
- [ ] **Nyquist-zeroing audit.** Construct a one-mode input that lives exactly at the x-Nyquist (`Nx/2`) and verify `spectral::partial<3, 0, 0>` returns the zero field (per the documented Nyquist-zeroing-for-odd-rank convention) — at $N = 16$ to keep the test fast.
- [ ] All Phase 1–8 tests still pass without modification.

## Documentation

- [ ] Every new public symbol carries a Doxygen block including `\brief`, `\param`, `\returns`, `\tparam` (where applicable), `\since 0.9.0`. Covers: `spectral::partial`, `spectral::laplacian`, `spectral::biharmonic`, `spectral::gradient`, `spectral::forwardR2C`, `spectral::backwardC2R`, `spectral::ComplexField`, `spectral::kxRfft`, `spectral::kyFull`, `spectral::kzFull`.
- [ ] Each spectral header carries the explicit "verification only — do not use as a production differentiation engine" admonition in its file-level Doxygen block, citing `mission.md`.
- [ ] User Guide note added at `docs/user-guide/notes/phase-9-spectral-verification.md` covering the hybrid API, verification-only positioning, build flag, and a worked CI-style assertion example.
- [ ] Developer Note added at `docs/developer-note/notes/phase-9-spectral-verification.md` with the five-section structure (Theory · Math derivation · Algorithm · Design decisions · References) per `specs/CLAUDE.md` step 4d. References section non-empty: at least one external citation (DOI / textbook bibliographic entry / permanent URL).
- [ ] `specs/tech-stack.md` updated: PocketFFT recorded as vendored under `extern/pocketfft/`; build flag documented.
- [ ] `CHANGELOG.md` entry under `0.9.0` lists the vendored PocketFFT, the new `gridcalc::spectral` namespace, the build flag, and the FD–FFT cross-check coverage.
- [ ] `STATUS.md` updated: Phase 9 row marked Done; Phase 10 stated as the new Next action. New "Decisions worth knowing" entries for vendored PocketFFT, hybrid spectral API, verification-only role, and Nyquist-zeroing convention.
- [ ] CI `render-docs` job produces fresh aggregate PDFs that include the Phase 9 chapter in both books.

## Performance

- N/A at Phase 9. The cross-check fixture runs ~70 tests at $N \in \{16, 32, 64, 128\}$. Estimated CI runtime: ~90–120 s on Ubuntu Clang. Phase 20 owns optimization; the Phase-9 spectral path is verification-only and not on a perf-sensitive code path.
