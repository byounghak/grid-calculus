# Plan: Phase 4 — Functional evaluation, callable API (scalar, periodic)

## Group 1 — Add `func::evaluate` with SFINAE-detected arity

- Create `include/gridcalc/func/Functional.hpp` (filename matches the roadmap; namespace stays `gridcalc::func`).
- Define `template <class F, class Tag = Pairwise> double evaluate(const Field<double>& psi, F&& f, Tag tag = {})` using `if constexpr` + `std::is_invocable_r_v<double, F&, ...>` to dispatch on the callable's arity:
  - `f(double)` — only `ψ` needed; skip `gradient` and `laplacian`.
  - `f(double, const Vec3d&)` — `ψ` and `∇ψ` needed; skip `laplacian`.
  - `f(double, const Vec3d&, double)` — full-arity; materialize both `gradient` and `laplacian`.
  - Otherwise: `static_assert` via a `detail::deferred_false<F>` trick to surface a clean error message at the call site.
- Implementation: eager materialization. Internally calls `diff::gradient(psi)` / `diff::laplacian(psi)` only as required by the detected arity; then loops over grid points filling a `Field<double>` integrand; then calls `func::integrate(integrand, tag)`.
- `Tag` is constrained to `Pairwise` or `Kahan` (the existing tags from Phase 3 `Integrate.hpp`); a `static_assert` on the tag type guards against accidental misuse.
- No source files added; `gridcalc` stays INTERFACE.
- Doxygen blocks on every public symbol with `\since 0.4.0`.

**Commit message:** `feat(func): add evaluate() with SFINAE-detected callable arity`

## Group 2 — Tests for evaluate()

- Add `test/functional_test.cpp`:
  - `FunctionalTest.GinzburgLandauHandComputed` — at `N = 64`, ψ = sin(x) sin(y) sin(z) on `[0, 2π]³`, W = 1. Hand-computed reference: `F = (3/2)π³ + (411/64)π³ = (507/64)π³ ≈ 245.6`. Pass when relative error < 1e-2 (gradient term dominates the discrete error and is `O(h²)`-accurate via Phase 1+2 stencils).
  - `FunctionalTest.GinzburgLandauSecondOrderConvergence` — sweep `N ∈ {16, 32, 64}`, log-log slope of relative error in `[1.8, 2.2]`, halving-h ratio in `[3.5, 4.5]`.
  - `FunctionalTest.PsiOnlyArity` — `f(ψ) = ψ²`, hand-computed reference π³ for ψ = sin(x) sin(y) sin(z) over `[0, 2π]³`. Verifies the 1-arg dispatch path skips both `gradient` and `laplacian`.
  - `FunctionalTest.PsiAndGradArity` — `f(ψ, ∇ψ) = (1/2)|∇ψ|²`, reference `(3/2)π³`. Verifies the 2-arg dispatch path skips only `laplacian`.
  - `FunctionalTest.PsiGradLapArity` — `f(ψ, ∇ψ, ∇²ψ) = ψ * ∇²ψ`, hand-computed reference `-3π³` for ψ = sin(x) sin(y) sin(z) (since ∇²ψ = -3ψ). Verifies the 3-arg dispatch path materializes both gradient and laplacian.
  - `FunctionalTest.KahanTagAgrees` — same GL setup, both reducer tags produce results agreeing within 1e-13 relative.
- Wire `functional_test.cpp` into `test/CMakeLists.txt`.

**Commit message:** `test(func): add evaluate accuracy + arity-dispatch tests`

## Group 3 — Phase 4 doc-notes (User Guide + Developer Note)

- `docs/user-guide/notes/phase-4-functional-evaluation.md` — user-facing summary: signature of `evaluate`, the three accepted callable arities and what each means semantically, the role of the reducer tag (forwarded to `integrate`), and a worked example that walks through the Ginzburg–Landau free energy from setup through hand-computed reference. This is the "Tutorial doc page walking through the example" deliverable from `roadmap.md` Phase 4.
- `docs/developer-note/notes/phase-4-functional-evaluation.md` — five-section structure per `specs/CLAUDE.md` step 4d:
  1. **Theory**: the variational object `F[ψ] = ∫ f(ψ, ∇ψ, ∇²ψ) dV`; periodic Cartesian setting; how the discrete `F_h` inherits its convergence from the underlying stencils (gradient is `O(h²)`, Laplacian is `O(h²)`, integration is spectrally accurate, so the bottleneck is whichever derivative the integrand actually uses).
  2. **Math derivation**: leading-order discrete error for the GL example; explicit closed-form for `F` on ψ = sin(x) sin(y) sin(z) including the `(411/64)π³` decomposition; argument that `evaluate` is exact-to-machine-precision when `f` does not involve any derivative.
  3. **Algorithm**: the SFINAE arity-dispatch ladder; eager materialization rationale; reuse of Phase 1–3 primitives (`diff::gradient`, `diff::laplacian`, `func::integrate`); algorithmic complexity.
  4. **Design decisions**: why SFINAE-detected arity (cite `requirements.md`); why eager rather than per-point lazy; why the worked example is GL.
  5. **References**: Provatas & Elder, *Phase-Field Methods in Materials Science and Engineering* (textbook citation, GL functional in chapter 3); Boyd, *Chebyshev and Fourier Spectral Methods*, 2nd ed. (textbook, on tensor-product convergence); the `std::is_invocable_r` cppreference page (permanent URL).

**Commit message:** `docs(notes): add Phase 4 User Guide and Developer Note notes`

## Group 4 — Bookkeeping

- Bump `CMakeLists.txt` `project(gridcalc VERSION 0.3.0)` → `0.4.0`.
- Update `test/version_test.cpp` (`VersionTest.ReturnsZeroThreeZero` → `ReturnsZeroFourZero`).
- `CHANGELOG.md`: new `## 0.4.0 — Phase 4: functional evaluation` entry.
- `specs/STATUS.md`: bump "Last updated"; move Phase 4 to Done; update "Where the project stands" + "Next action" (which becomes Phase 5 — Explicit Euler diffusion solver).

**Commit message:** `chore: bump version to 0.4.0 and refresh STATUS for Phase 4`
