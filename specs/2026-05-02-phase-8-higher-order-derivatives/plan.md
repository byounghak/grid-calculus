# Plan: Phase 8 — Higher-order spatial derivatives (3rd, 4th)

## Group 1 — Add `stencil::ThirdDerivative<Order>` and `stencil::FourthDerivative<Order>`

- Create `include/gridcalc/stencil/ThirdDerivative.hpp` mirroring `FirstDerivative.hpp`'s pattern:
  - Primary template `template <int Order> struct ThirdDerivative;` left undefined.
  - `ThirdDerivative<2>`: `radius = 2`, `offsets = {-2, -1, 0, 1, 2}`, `weights = {-1.0/2, 1.0, 0.0, -1.0, 1.0/2}`. Unscaled — consumer multiplies by `1/h^3`. (Anti-symmetric stencil; odd derivative.)
  - `ThirdDerivative<4>`: `radius = 3`, `offsets = {-3, -2, -1, 0, 1, 2, 3}`, `weights = {1.0/8, -1.0, 13.0/8, 0.0, -13.0/8, 1.0, -1.0/8}`. Unscaled — consumer multiplies by `1/h^3`.
  - Doxygen blocks include `\since 0.8.0`.
- Create `include/gridcalc/stencil/FourthDerivative.hpp` mirroring `CentralDifference.hpp`'s pattern:
  - Primary template undefined.
  - `FourthDerivative<2>`: `radius = 2`, `offsets = {-2, -1, 0, 1, 2}`, `weights = {1.0, -4.0, 6.0, -4.0, 1.0}`. Unscaled — consumer multiplies by `1/h^4`. (Symmetric; even derivative.)
  - `FourthDerivative<4>`: `radius = 3`, `offsets = {-3, -2, -1, 0, 1, 2, 3}`, `weights = {-1.0/6, 2.0, -13.0/2, 28.0/3, -13.0/2, 2.0, -1.0/6}`. Unscaled — consumer multiplies by `1/h^4`.
  - Doxygen blocks include `\since 0.8.0`.
- Both files follow the existing code style (`\since`, `\tparam`, unitless integer offsets, `static constexpr` arrays).
- The primary templates remain undefined; instantiating `ThirdDerivative<3>` / `FourthDerivative<6>` etc. continues to produce a compile error at the call site.

**Commit message:** `feat(stencil): add ThirdDerivative<Order> and FourthDerivative<Order> tables`

## Group 2 — Implementation helper for tensor-product partial derivatives

- Create `include/gridcalc/diff/detail/MultiIndexDerivative.hpp` (a `detail::` namespace under `gridcalc::diff`):
  - `template <int Nx, int Ny, int Nz, int Order> core::Field<double> mixedDerivative(const core::Field<double>& psi)` — computes $\partial^{Nx+Ny+Nz}\psi / \partial x^{Nx}\partial y^{Ny}\partial z^{Nz}$ in a single pass over the grid via the outer product of three 1D weight tables.
  - Picks the per-axis 1D weight table by `constexpr if`-equivalent (template specialization or `if constexpr`):
    - `Na = 0`: identity (single-element table `{1.0}`, offset `{0}`, no scaling).
    - `Na = 1`: `stencil::FirstDerivative<Order>`, scaled by `1/h_a`.
    - `Na = 2`: `stencil::Coefficients<Order>`, scaled by `1/h_a^2`.
    - `Na = 3`: `stencil::ThirdDerivative<Order>`, scaled by `1/h_a^3`.
    - `Na = 4`: `stencil::FourthDerivative<Order>`, scaled by `1/h_a^4`.
  - Outer-product loop nest: triple-nested over `(mx, my, mz)` ∈ axis-x-table × axis-y-table × axis-z-table at each grid point `(i, j, k)`. The outer-product weight is the product of the three 1D weights; the read is `field(i+offset_x, j+offset_y, k+offset_z)`. Periodic wrap delegated to `Field::Policy`.
  - Function is `inline` and lives in a header (consistent with the rest of `diff/`).
- This helper is the single implementation site for all 31 named partial-derivative functions in Groups 3, 4, 5. The named functions are 2-line wrappers.

**Commit message:** `feat(diff): add multi-index tensor-product derivative helper`

## Group 3 — Rank-2 partials in `diff/MixedPartial.hpp`

- Create `include/gridcalc/diff/MixedPartial.hpp` with six functions, each templated on `<int Order = 2>`:
  - `dxx<Order>(psi)` → `detail::mixedDerivative<2, 0, 0, Order>(psi)`
  - `dyy<Order>(psi)` → `detail::mixedDerivative<0, 2, 0, Order>(psi)`
  - `dzz<Order>(psi)` → `detail::mixedDerivative<0, 0, 2, Order>(psi)`
  - `dxy<Order>(psi)` → `detail::mixedDerivative<1, 1, 0, Order>(psi)`
  - `dxz<Order>(psi)` → `detail::mixedDerivative<1, 0, 1, Order>(psi)`
  - `dyz<Order>(psi)` → `detail::mixedDerivative<0, 1, 1, Order>(psi)`
- Each carries a Doxygen block with `\brief`, `\param`, `\returns`, `\tparam Order`, `\since 0.8.0`.
- The header file name (`MixedPartial.hpp`) matches the roadmap deliverable text; the *symbol* names follow the project's d-prefix convention.

**Commit message:** `feat(diff): add rank-2 partial derivatives (dxx, dyy, dzz, dxy, dxz, dyz)`

## Group 4 — Rank-3 partials in `diff/ThirdOrder.hpp`

- Create `include/gridcalc/diff/ThirdOrder.hpp` with ten functions, each templated on `<int Order = 2>`:
  - Single-axis: `d3dx3`, `d3dy3`, `d3dz3`.
  - Two distinct axes (multiplicity (2, 1)): `d3dx2dy`, `d3dxdy2`, `d3dx2dz`, `d3dxdz2`, `d3dy2dz`, `d3dydz2`.
  - Three distinct axes: `d3dxdydz`.
- Each is a 2-line wrapper around `detail::mixedDerivative<Nx, Ny, Nz, Order>`.
- Doxygen blocks include `\since 0.8.0`.

**Commit message:** `feat(diff): add rank-3 partial derivatives`

## Group 5 — Rank-4 partials in `diff/FourthOrder.hpp` and biharmonic in `diff/Biharmonic.hpp`

- Create `include/gridcalc/diff/FourthOrder.hpp` with fifteen functions:
  - Single-axis: `d4dx4`, `d4dy4`, `d4dz4`.
  - Two distinct axes, multiplicity (3, 1): `d4dx3dy`, `d4dxdy3`, `d4dx3dz`, `d4dxdz3`, `d4dy3dz`, `d4dydz3`.
  - Two distinct axes, multiplicity (2, 2): `d4dx2dy2`, `d4dx2dz2`, `d4dy2dz2`.
  - Three distinct axes, multiplicity (2, 1, 1): `d4dx2dydz`, `d4dxdy2dz`, `d4dxdydz2`.
- Each is a 2-line wrapper around `detail::mixedDerivative<Nx, Ny, Nz, Order>`. `\since 0.8.0`.
- Create `include/gridcalc/diff/Biharmonic.hpp`:
  - `template <int Order = 2> inline core::Field<double> biharmonic(const core::Field<double>& psi)` — sums six tensor-product partial derivatives in a single fused pass over the grid (not via `laplacian(laplacian(psi))`):
    $\nabla^4 \psi = \partial_x^4 \psi + \partial_y^4 \psi + \partial_z^4 \psi + 2(\partial_x^2 \partial_y^2 \psi + \partial_x^2 \partial_z^2 \psi + \partial_y^2 \partial_z^2 \psi)$.
  - Implementation: extend the `detail::mixedDerivative` helper pattern to accumulate six weighted contributions into one result `Field` in a single pass — or, equivalently, call `detail::mixedDerivative` six times and add. Pick whichever keeps the loop nest readable; the perf difference at Phase 8's stage is negligible (Phase 20 owns optimization).
  - `template <int Order = 2> inline core::Field<double> d4(const core::Field<double>& psi) { return biharmonic<Order>(psi); }` — the d-prefix-family alias.
  - Doxygen blocks include `\since 0.8.0`. The alias's block notes that it forwards to `biharmonic<Order>`.

**Commit message:** `feat(diff): add rank-4 partials and biharmonic with d4 alias`

## Group 6 — Tests

- Add `test/phase_8_test.cpp` containing:
  - **Weight-table audits.** `EXPECT_DOUBLE_EQ` on every weight in `ThirdDerivative<2>`, `ThirdDerivative<4>`, `FourthDerivative<2>`, `FourthDerivative<4>` against the hard-coded reference values from the Developer Note's Math derivation section.
  - **Convergence-sweep fixture.** A typed/parameterized GoogleTest fixture `HigherOrderDerivativeTest` that, for each (operator, Order) pair, runs the sweep `N ∈ {16, 32, 64}` on the periodic `[0, 2π)^3` box with $\psi = \sin(x)\sin(y)\sin(z)$, computes max-norm error against the analytical partial derivative (closed-form trig), and asserts the log-log slope lies in `[Order - 0.5, Order + 0.5]`. Operator coverage: all 31 named partials at `Order ∈ {2, 4}` = 62 sweeps, plus `biharmonic<2>` and `biharmonic<4>` = 64 sweeps total.
    - Implementation: a tag struct per operator (`struct OpDxxOrder2 { static constexpr int Nx=2, Ny=0, Nz=0, Order=2; static auto apply(const auto& f) { return diff::dxx<2>(f); } };` etc.), wired into `TYPED_TEST_SUITE`.
    - Saves writing 62 separate `TEST(...)` macros while keeping each (op, Order) pair as its own observable test result.
  - **Roadmap acceptance — biharmonic recovers $|\mathbf{k}|^4$.** `BiharmonicTest.RecoversKToTheFourthOrder2` and `BiharmonicTest.RecoversKToTheFourthOrder4`: for $\psi = \sin(\mathbf{k}\cdot\mathbf{x})$ on `[0, 2π)^3`, two non-trivial integer wave vectors (e.g., $\mathbf{k} = (1, 1, 1)$ giving $|\mathbf{k}|^4 = 9$ and $\mathbf{k} = (1, 2, 1)$ giving $|\mathbf{k}|^4 = 36$), `biharmonic<Order>(psi)` recovers $|\mathbf{k}|^4 \psi$ pointwise. Tolerance is the order-`Order` truncation error at `N = 64`.
- Wire `test/phase_8_test.cpp` into `test/CMakeLists.txt` alongside the existing `*_test.cpp` files.
- Existing tests (Phase 1–7) are untouched.

**Commit message:** `test(diff): convergence-order sweeps for all Phase 8 partials and biharmonic acceptance`

## Group 7 — Phase 8 doc-notes

- `docs/user-guide/notes/phase-8-higher-order-derivatives.md` — user-facing summary:
  - The full list of 31 named partial-derivative functions plus `biharmonic` / `d4`.
  - Naming convention rules in one paragraph (rank prefix + lexicographic axis order + omit exponent 1).
  - Calling convention (`auto h = diff::d3dx2dy<4>(psi);` etc.).
  - The order-template parameter mirrors Phase 7.
  - One worked example: a 3D Laplacian-of-Laplacian regularizer ($\nabla^4$) on a Gaussian-bump field, computed via `biharmonic<4>(psi)` and shown to agree with $|\mathbf{k}|^4$ in spectral expectation.
- `docs/developer-note/notes/phase-8-higher-order-derivatives.md` — five-section structure per `specs/CLAUDE.md` step 4d:
  1. **Theory.** The space of multi-indices $\boldsymbol\alpha = (n_x, n_y, n_z)$ with $|\boldsymbol\alpha| \leq 4$ on a Cartesian-orthogonal periodic grid; tensor-product structure of the partial-derivative operator on a separable grid; how it lets the multi-axis stencil factor into outer products of 1D stencils. The biharmonic as a contracted scalar operator $\sum_{i,j} \partial_i^2 \partial_j^2$.
  2. **Math derivation.** Explicit Taylor-matching derivation of the four new 1D weight tables:
     - `ThirdDerivative<2>::weights = {-1/2, 1, 0, -1, 1/2}` (5-point anti-symmetric).
     - `ThirdDerivative<4>::weights = {1/8, -1, 13/8, 0, -13/8, 1, -1/8}` (7-point anti-symmetric).
     - `FourthDerivative<2>::weights = {1, -4, 6, -4, 1}` (5-point symmetric).
     - `FourthDerivative<4>::weights = {-1/6, 2, -13/2, 28/3, -13/2, 2, -1/6}` (7-point symmetric).
     Show that an outer product of 1D weight tables yields the correct multi-axis partial weights to the *minimum* of the per-axis truncation orders. Closed-form analytical reference for $\partial^{\boldsymbol\alpha} \psi$ on $\psi = \sin(x)\sin(y)\sin(z)$, used as the manufactured-solution reference for the convergence sweep.
  3. **Algorithm.** The `detail::mixedDerivative<Nx, Ny, Nz, Order>` helper: which 1D table is selected per axis, how the outer-product loop nest scales (per-grid-point work is $\prod_a (2 r_a + 1)$), file-path references to `include/gridcalc/diff/detail/MultiIndexDerivative.hpp`, the named-function wrappers, and the biharmonic's six-term fused implementation. Tests live in `test/phase_8_test.cpp`.
  4. **Design decisions.** One short paragraph each, citing `requirements.md`:
     - Scope expansion to all 31 unique partials.
     - d-prefix naming convention with lexicographic axis order.
     - Direct stencils over `laplacian∘laplacian` composition for the biharmonic (cite Phase 2's `divergence(gradient)` stride-2 issue).
     - Sibling templates over a unified `Derivative<Rank, Order>` redesign.
     - Hard-coded weights over a Fornberg generator (cite Phase 7's same call).
  5. **References.** At least: Fornberg (1988) for the weight-derivation algorithm (DOI). LeVeque (2007) for the standard finite-difference table (full bibliographic entry). NIST DLMF or a stable Wikipedia URL for the biharmonic operator definition. Cross-references to Phase 7's Developer Note count as internal context but do not satisfy the "external citation" minimum.
- Both notes follow the structure documented in `docs/user-guide/notes/README.md` and `docs/developer-note/notes/README.md`.

**Commit message:** `docs(notes): add Phase 8 User Guide and Developer Note notes`

## Group 8 — Bookkeeping

- Bump `CMakeLists.txt` `project(gridcalc VERSION 0.7.0)` → `0.8.0`.
- Update `test/version_test.cpp` (`ReturnsZeroSevenZero` → `ReturnsZeroEightZero`).
- `CHANGELOG.md`: new `## 0.8.0 — Phase 8: Higher-order spatial derivatives` entry listing the new stencil tables, the 31 named partial-derivative functions, `biharmonic` + `d4` alias, and the convergence-test coverage.
- `specs/STATUS.md`: bump "Last updated"; move Phase 8 row to Done; update "Where the project stands" + "Next action" (Phase 9 — spectral verification harness). Add new "Decisions worth knowing" entries for the d-prefix naming convention and for the second deferral of the unified `Derivative<Rank, Order>` template.

**Commit message:** `chore: bump version to 0.8.0 and refresh STATUS for Phase 8`
