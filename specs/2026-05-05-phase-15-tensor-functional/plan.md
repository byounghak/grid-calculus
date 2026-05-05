# Plan: Phase 15 — Functional evaluation for vector/tensor data

Six commit-sized groups. Build and tests stay green between every commit.

## Group 1 — Spec files

- `specs/2026-05-05-phase-15-tensor-functional/{requirements, plan, validation}.md`.

**Commit message:** `docs: phase 15 spec files (tensor functional + ET layer)`

## Group 2 — Contraction expression-template scaffolding

- New directory `include/gridcalc/tensor/expr/` and umbrella header `include/gridcalc/tensor/Expressions.hpp`. Per-node headers under `tensor/expr/`: `Leaf.hpp`, `IdentityField.hpp`, `Scaled.hpp`, `Plus.hpp`, `Minus.hpp`, `Negate.hpp`, `Trace.hpp`, `Sym.hpp`, `SingleContract.hpp`, `DoubleContract.hpp`. Each node exposes:
  - `using value_type = ...;` — `double`, `core::Vec3d`, or `core::Mat3d`.
  - `const core::Grid& grid() const noexcept;` — propagated up the tree.
  - `value_type evalAt(int i, int j, int k) const noexcept;` — the per-point evaluator that the materialize / reduce passes call.
- Operator overloads in the umbrella header: unary `-`, binary `+` / `-` (same-`value_type`), scalar `*` / `/` (commutative, both sides). A `Field<T>` participates implicitly via overload resolution that wraps it in `expr::Leaf<T>`.
- Phase 13's `core::Field<T>` does not need to change; the ET nodes hold a `const Field<T>&` reference.
- `tensor::expr::field(const Field<T>&)`, `tensor::expr::identityField(const Grid&)` factory helpers.
- `tensor::expr::materialize(E)` (works for any expression node) returns `core::Field<value_type<E>>`. Allocates one Field, single-pass loop over all cells. Used as an explicit escape hatch and as the implementation of the implicit `operator core::Field<value_type<E>>()` conversion on each node.
- `func::integrate(E expr, Tag tag = Pairwise{})` overload added to `include/gridcalc/func/Integrate.hpp` (or a new sibling header `include/gridcalc/func/IntegrateExpr.hpp` if separating keeps Phase 3's existing surface clean — choice made during implementation, with a preference for adding to the existing `Integrate.hpp` to minimize new surface area). Requires `value_type<E>` to be `double`. Walks all cells once, evaluates `expr.evalAt(i, j, k)`, accumulates with the existing `Pairwise` / `Kahan` reduction. Bit-identical to `func::integrate(materialize(expr), tag)` save the temporary allocation.
- New file `test/tensor_expressions_test.cpp` covering:
  - **Materialization equivalence**: each ET node's `materialize(expr)` matches Phase 13's eager `tensor::trace`/`singleContract`/`doubleContract` byte-for-byte (max-abs `< 1e-15`) on a manufactured trig field.
  - **Operator algebra**: `(a + b)` materializes to the pointwise sum of `materialize(a)` and `materialize(b)`; same for `-`, `scalar *`, `scalar /`.
  - **Identity-field correctness**: `identityField(grid).evalAt(i, j, k) == Mat3d::Identity()` at every cell.
  - **`sym(M)` correctness**: matches `0.5 * (M + Mᵀ)` on a non-symmetric input.
  - **Fused integrate**: `func::integrate(0.5 * doubleContract(a, b))` matches `0.5 * sum_over_cells(doubleContract_eager(a, b)) * dV` with no allocation of an intermediate `Field<double>` (verified via assertion that `materialize` is *not* called on the doubleContract node — checked via a counter on a debug-only ET node, or simply by an algorithmic-cost test).
  - **Mixed-type compositions**: `(scalar * trace(M)) + (scalar * doubleContract(A, B))` evaluates correctly across nested levels.

  Tests target ~14–20 cases.
- Doxygen blocks on every new public symbol with `\since 0.15.0`. Internal nodes under `tensor/expr/` are still public-facing (per the project's "don't put it under `detail/` unless truly internal" rule), so they get full Doxygen too.

**Commit message:** `feat(tensor): contraction expression-template layer (lazy nodes + fused integrate)`

## Group 3 — `func::evaluate` overloads for `Field<Vec3d>`

- Extend `include/gridcalc/func/Functional.hpp` with two new SFINAE-detected overloads:
  - `template <class F, class Tag = Pairwise> double evaluate(const Field<Vec3d>& u, F&& f, Tag tag = {})`
  - Internal dispatch ladder: longest-arity wins.
    - `f(const Vec3d& u, const Mat3d& grad_u)` → materialize `diff::gradient(u)` once, loop pointwise.
    - `f(const Vec3d& u)` → no derivatives; loop pointwise.
- `static_assert` failure path in the same shape as the scalar overload.
- Doxygen block on the new overload calling out the rank-2 gradient convention `M(i, j) = ∂_j u_i` and `\since 0.15.0`.
- New tests in a new file `test/func_evaluate_vector_test.cpp`:
  - `EvaluateUOnlyArity` — `f(u) = ‖u‖²` integrated against `u = (sin x, cos y, sin z cos x)` matches the closed form.
  - `EvaluateUAndGradArity` — `f(u, ∇u) = trace(∇u)²` integrated against `u = (sin x, sin y, sin z)` matches `∫ (cos x + cos y + cos z)² dV` to `O(h²)` convergence.
  - `LongestArityWins` — passing a 2-arg-callable that's also implicitly callable as 1-arg picks the 2-arg branch.
  - `RejectsBadCallable` — passing a callable with the wrong signature (e.g., taking `Mat3d`) triggers the `static_assert`. Verified by `static_assert`-presence comment + a note that the failure is compile-time-only (no runtime test).
  - `KahanReductionTagWorks` — same integrand, `Pairwise` and `Kahan` tags agree to round-off.
  - ~6 tests total.

**Commit message:** `feat(func): evaluate over Field<Vec3d> with (u) and (u, grad_u) arities`

## Group 4 — Elastic-energy shared helpers + acceptance test

- New header `examples/common/elastic_energy.hpp` exposing:
  - `Field<Vec3d> makeUniaxialPeriodicDisplacement(const Grid& grid, double alpha, int axis = 0)` — periodic uniaxial dilation wave `u_axis(x_axis) = α sin(k₀ x_axis)` with `k₀ = 1` (the box is `[0, 2π]³` per the project convention; `axis ∈ {0, 1, 2}` selects which Cartesian direction the dilation is along).
  - `inline double linearElasticEnergyDensity(double lambda, double mu, const Mat3d& grad_u)` — pointwise `½ σ:ε` with `ε = sym(grad_u)`, `σ = λ tr(ε) I + 2μ ε`, returning `½ (λ tr(ε)² + 2μ ε:ε)`. The single-line Eigen expression for performance and clarity.
  - `inline double analyticalLinearElasticEnergy(const Grid& grid, double alpha, double lambda, double mu, int axis = 0)` — closed-form `2 π³ (λ + 2μ) α²` (independent of `axis` by symmetry; `axis` retained as a parameter for explicit symmetry).
  - `lameFromYoungPoisson(double E, double nu) -> std::pair<double, double>` returning `{lambda, mu}` for callers that prefer the engineering-constant inputs.
- New file `test/elastic_energy_test.cpp` discharging the roadmap acceptance bar:
  - `MatchesAnalyticalAtN64` — at `N = 64`, `α = 0.01`, `E = 1`, `ν = 0.3`: relative error `|F_h - F_ref| / |F_ref| < 1e-3`. Uses `func::evaluate(u, ...)` with the energy-density callable.
  - `SecondOrderConvergence` — sweep `N ∈ {16, 32, 64}`; ratio `e(N=16) / e(N=32)` and `e(N=32) / e(N=64)` both ≈ 4 (slope ∈ `[1.6, 2.4]` on `log(e) / log(h)`).
  - `EnergyDensityFormulaAtZeroStrain` — `linearElasticEnergyDensity(lambda, mu, Mat3d::Zero()) == 0`.
  - `EnergyDensityFormulaAtUnitStrainAlongX` — `linearElasticEnergyDensity(λ, μ, diag(1, 0, 0))` equals `½ (λ + 2μ)` analytically.
  - `EachAxisGivesSameAnalyticalEnergy` — uniaxial wave along `x`, `y`, `z` axes all give numerically equal energies (symmetry check).
  - `MatchesViaExpressionLayer` — alternative computation via the ET layer (`func::integrate(0.5 * (lambda * trace_expr(eps_expr)^2 + 2.0 * mu * doubleContract_expr(eps_expr, eps_expr)))`) matches the `func::evaluate` route to round-off.
  - ~6 tests.

**Commit message:** `feat(examples): elastic-energy helpers + acceptance test (rank-2 functional surface)`

## Group 5 — Demo binary + CI gate

- New file `examples/elastic_energy.cpp` — CLI demo. Options: `--n-x` (default 64), `--lambda` (default `λ(E=1, ν=0.3)`), `--mu` (default `μ(E=1, ν=0.3)`), `--alpha` (default 0.01), `--axis {x|y|z}` (default `x`), `--out-dir` (default `./out/elastic-energy/`). Computes the displacement, the strain tensor, the energy density, and the total energy; writes a single VTK snapshot `elastic_energy.vtk` containing the `ε_xx` field (or the diagonal axial-strain component for the chosen axis); prints to stdout the computed `F`, the analytical `F_ref`, the relative error, and the `λ`, `μ`, `α` parameters used. Wired into `examples/CMakeLists.txt`.
- New file `test/elastic_energy_demo_test.cpp` — small CI gate (~5 s budget) that runs a stripped-down version of the demo at `N = 24` and asserts the relative error is `< 1e-2` (loose since `N = 24` is below the convergence-test grid). Mirrors the Phase 12 / 14 static-cached pattern; no actual VTK output written from the test.
- New file `scripts/render_elastic_snapshots.py` (matplotlib `Agg` backend, mirroring `scripts/render_diffusion_snapshots.py`) renders three PNG snapshots from `α ∈ {0.005, 0.01, 0.02}` runs on `N = 64`, showing the `cos(k₀ x)` strain profile across the `(y = 0, z = 0)` slice. Snapshots committed under `docs/user-guide/figures/elastic-energy/{small, medium, large}.png`.

**Commit message:** `feat(examples): elastic-energy demo binary + CI gate`

## Group 6 — Doc deliverables and bookkeeping

- **Chapter slot bookkeeping.** `docs/user-guide/chapters/15-diamond-lattice.tex` is currently the placeholder for Phase 17 (it survived the 2026-05-04 renumber unchanged). Rename it to `docs/user-guide/chapters/17-diamond-lattice.tex` to free slot 15 for Phase 15. Update the `\input{}` line in `docs/user-guide/main.tex` accordingly. (Slot 16 stays empty, reserved for Phase 16's non-orthogonal-Bravais chapter.) Phase 14's previously-touched `docs/user-guide/chapters/15-diamond-lattice.tex` "Phase 16 → Phase 17" comment edit carries over; no further content edits to the diamond placeholder.
- New file `docs/user-guide/chapters/15-tensor-functional.tex` — full chapter 15. Sections: motivation (functionals over vector / tensor fields, the rank-2 case); `func::evaluate(Field<Vec3d>, ...)` API with the rank-2 gradient convention; the contraction ET layer (operator algebra, `materialize`, `func::integrate(expr)` fused reduction); the elastic-energy worked example with the periodic-uniaxial-wave manufactured solution and the closed-form derivation; demo walkthrough (Group 5) with the three committed PNG snapshots; what this does *not* do (rank-3 / rank-4 deferred, periodicity-only, etc.).
- New file `docs/developer-note/chapters/14-tensor-functional.tex` — five-section structure (Theory · Math derivation · Algorithm · Design decisions · References). Theory: tensor-rank arithmetic, Einstein-summation contractions, isotropic linear elasticity (Lamé form). Math derivation: $\boldsymbol{\sigma} = \lambda \mathrm{tr}(\boldsymbol{\varepsilon})\mathbf{I} + 2\mu \boldsymbol{\varepsilon}$ from minor-symmetric stiffness; closed-form energy of the periodic uniaxial wave; `O(h²)` convergence argument inheriting from `diff::gradient`. Algorithm: file layout, the ET node design, `evalAt` recursion, the fused `func::integrate(expr)` loop. Design decisions: the four AskUserQuestion answers + the periodic-uniaxial-wave substitution + the Phase-13-eager-API-stays-untouched decision. References: Gurtin 1981 (linear elasticity), Holzapfel 2000 (continuum mechanics), Vandevoorde / Josuttis / Gregor *C++ Templates* 2nd ed. ch. 27 (expression templates), Veldhuizen 1995 (DOI:10.1145/216973.216979) for the ET technique, Eigen project docs.
- New three PNG snapshots under `docs/user-guide/figures/elastic-energy/{small, medium, large}.png`.
- Update `docs/user-guide/main.tex` and `docs/developer-note/main.tex` with `\input{}` lines for the new chapters.
- Bump `project(... VERSION 0.14.5)` → `0.15.0` in root `CMakeLists.txt`.
- Update `test/version_test.cpp` to assert `0.15.0`.
- New `## 0.15.0 — Phase 15` block in `CHANGELOG.md`.
- Refresh `specs/STATUS.md`: Phase 15 row → Done; "Last updated" date bump; "Next action" rewritten for Phase 16 (Non-orthogonal Bravais lattices); new decisions added.

**Commit message:** `chore: bump version to 0.15.0 and refresh STATUS for Phase 15`

---

After Group 6: `git push -u origin 2026-05-05-phase-15-tensor-functional`, open PR, watch CI per `specs/CLAUDE.md` step 6.
