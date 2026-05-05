# Requirements: Phase 15 — Functional evaluation for vector/tensor data

## Roadmap reference

> **Goal.** Functionals over vector and tensor fields.
>
> **Deliverables.**
>
> - Generalized `Functional` that accepts callables over `(T, Gradient<T>, ...)` argument packs.
> - Worked example: linear elastic energy $F = \int \tfrac{1}{2}\,\boldsymbol{\sigma}:\boldsymbol{\varepsilon}\,dV$.
>
> **Acceptance.** Elastic energy of a simple stretched cubic block matches the analytical value.

(Phase 15 in `specs/roadmap.md` after the 2026-05-04 Phase 14 (FVM) renumber; was the original Phase 14 before that insertion.)

## Decisions made for this feature

- **`func::evaluate` gains overloads on `Field<core::Vec3d>` input** with two SFINAE-detected callable arities:
  - `f(const core::Vec3d& u)` — only `u` is needed; no derivative is computed.
  - `f(const core::Vec3d& u, const core::Mat3d& grad_u)` — `u` and its rank-2 gradient `M(i, j) = ∂_j u_i` (Phase 13 convention; `trace = divergence` by construction).

  The longest-arity-wins tie-break and the `static_assert` failure path mirror the existing scalar `evaluate` (Phase 4 / Phase 11) signature exactly. Return type must be convertible to `double`. Reduction tag (`Pairwise` / `Kahan`) defaults to `Pairwise{}` per Phase 3.

- **Rank-3 / rank-4 tensor types are deferred.** The elastic-energy worked example with `ε = sym(∇u)` only needs `Field<core::Vec3d>` input and `Field<core::Mat3d>` derivatives — both already shipped in Phase 13. No rank-3 (`gradient(Field<Mat3d>)`) or rank-4 (constitutive stiffness `C_ijkl`) types are added in Phase 15. The `Eigen::TensorFixedSize` portability risk recorded under STATUS.md "Open / deferred items" stays deferred — a phase that genuinely consumes rank-3 data will pin the choice.

  Consequence: `func::evaluate` does **not** gain a `Field<core::Mat3d>` input overload in Phase 15, because the natural callable signature would take `(const Mat3d& m, const Tensor3& grad_m)` and rank-3 doesn't exist yet. Adding it later is additive (new overload, no signature change to the Phase 15 surface).

- **Contraction expression-template (ET) layer lands in Phase 15** (deferred from Phase 11 and Phase 13 per the standing "Decisions worth knowing" line). New header `include/gridcalc/tensor/Expressions.hpp` plus an internal `tensor/expr/` directory hosting the per-operator node types. Scope:

  - **Leaf / wrapper nodes.** `tensor::expr::field(const core::Field<T>&) -> Leaf<T>` for `T ∈ {double, core::Vec3d, core::Mat3d}`. `tensor::expr::identityField(const core::Grid&) -> IdentityField` broadcasts `Mat3d::Identity()` at every cell with no storage allocation.
  - **Pointwise arithmetic.** Operator overloads `+`, `-`, unary `-`, scalar `*` and `/` produce expression nodes. Same value type on both operands of `+` / `-`; mixed scalar ops broadcast. Leaf nodes participate in the algebra without an explicit `field(...)` wrap (a `Field<T>` is implicitly convertible to its leaf wrapper via overload resolution on the operator).
  - **Tensor primitives, lazy.** `tensor::expr::trace(E) -> double-valued expression` (where `E` produces `Mat3d`); `tensor::expr::sym(E) -> Mat3d-valued expression` returning `½(M + Mᵀ)`; `tensor::expr::singleContract(L, R) -> Mat3d-valued expression` (matrix product); `tensor::expr::doubleContract(L, R) -> double-valued expression` (full contraction `A : B = A_ij B_ij`).
  - **Materialization on demand.** `tensor::expr::materialize(E) -> Field<value_type<E>>` walks all grid points once, evaluates the expression's `evalAt(i, j, k)`, stores into a fresh `Field`. Implicit conversion from any expression node to `Field<value_type<E>>` is provided so existing Phase 13 call sites (`Field<double> tr = tensor::trace(M);`) continue to compile if they migrate from the eager `tensor::*` API to the lazy `tensor::expr::*` API.
  - **Fused reduction.** `func::integrate(expr, tag = Pairwise{})` overload (in `include/gridcalc/func/Integrate.hpp` or a new `tensor/Reduce.hpp`) consumes a double-valued expression and reduces in a single fused pass with no intermediate `Field` allocation. This is the primary motivation for the ET layer and the path the elastic-energy demo uses.

  **Scope boundary.** ET layer is fully-fused for the elastic-energy demo's Mat3d arithmetic: `½ (λ tr(ε)² + 2μ ε:ε)` evaluates with zero `Field<Mat3d>` temporaries. No fancy alias detection (the demo never aliases inputs); no SIMD / OpenMP retrofitting — those stay in Phase 21's performance pass.

- **Phase 13's eager `tensor::*` primitives stay untouched.** `tensor::trace`, `tensor::singleContract`, `tensor::doubleContract` continue to return materialized `Field<...>` (byte-for-byte unchanged). The new lazy versions live under `tensor::expr::*` and are a parallel surface, not a replacement. Phase 13's tests keep passing without modification.

- **Elastic-energy worked example ships as both a CI-gated acceptance test and a runnable demo binary** (mirroring Phase 12 CH and Phase 14 heterogeneous-diffusion). Shared helper `examples/common/elastic_energy.hpp` exposes:
  - `makeUniaxialPeriodicDisplacement(grid, alpha, axis = 0) -> Field<Vec3d>` — periodic uniaxial dilation wave (definition below).
  - `linearElasticEnergyDensity(lambda, mu, grad_u) -> double` — pointwise `½ σ:ε` returning a scalar, suitable as the `func::evaluate` callable on a `Field<Vec3d>` input.
  - `analyticalLinearElasticEnergy(grid, alpha, lambda, mu, axis) -> double` — closed-form reference for the periodic uniaxial wave.

- **Demo manufactured solution is a periodic uniaxial dilation wave.** Closed-form derivation:

  $$\mathbf{u}(x, y, z) = \bigl(A \sin(k_0\,x),\ 0,\ 0\bigr),\quad A = \alpha,\ k_0 = 1$$

  on $[0, 2\pi]^3$ (the project's standard periodic box). Then $\partial_x u_x = A\cos(k_0 x)$ and every other $\partial_j u_i = 0$, so the velocity-gradient $M(i, j) = \partial_j u_i$ is `diag(A cos(k₀ x), 0, 0)`. Strain $\boldsymbol{\varepsilon} = \mathrm{sym}(\nabla \mathbf{u}) = \nabla \mathbf{u}$ (already symmetric); stress $\boldsymbol{\sigma} = \lambda \mathrm{tr}(\boldsymbol{\varepsilon})\,\mathbf{I} + 2\mu\,\boldsymbol{\varepsilon}$. Energy density $\tfrac12 \boldsymbol{\sigma}:\boldsymbol{\varepsilon} = \tfrac12 (\lambda + 2\mu)\,A^2 \cos^2(k_0 x)$. Closed-form total:

  $$F_{\mathrm{ref}} = \tfrac12 (\lambda + 2\mu)\,A^2 \int_0^{2\pi}\!\cos^2(k_0 x)\,dx \cdot (2\pi)^2 = 2\pi^3\,(\lambda + 2\mu)\,A^2.$$

  **Why a sin-wave displacement and not a literal `u_x = α x` rigid stretch.** The roadmap acceptance bar says "a simple stretched cubic block." A literal rigid uniaxial stretch `u_x = α x` is **not** periodic on `[0, 2π]³` — it jumps from `2πα` to `0` across the periodic boundary, breaking the only `IndexPolicy::Periodic` shipped through Phase 14. The periodic uniaxial dilation wave above is the natural periodic analog: it has a constant *peak* strain `α` modulated sinusoidally across the box, and its closed-form energy admits the same $(\lambda + 2\mu)$ tensor signature as the literal rigid stretch (the trace and double-contract identities are exercised identically). Phase 19's non-periodic boundary-policy work will unlock the literal rigid-stretch geometry; until then the periodic wave is the right substitute. This decision is documented in the User Guide chapter 15 and the Developer Note chapter 14 so future readers don't re-litigate it.

- **Demo material parameters.** Isotropic linear elastic with Young's modulus `E = 1` and Poisson's ratio `ν = 0.3`. Lamé constants:
  - $\lambda = E\,\nu / ((1+\nu)(1-2\nu)) \approx 0.576923$
  - $\mu = E / (2(1+\nu)) \approx 0.384615$

  Demo amplitude `α = A = 0.01` (1 % peak strain — clearly within linear-elasticity validity). Closed-form total $F_{\mathrm{ref}} \approx 2\pi^3 \cdot 1.3461 \cdot 10^{-4} \approx 8.349 \times 10^{-3}$. Convergence sweep on $N \in \{16, 32, 64\}$.

- **Convergence acceptance bar.** `O(h²)` slope ∈ `[1.6, 2.4]` on the energy error $|F_h - F_{\mathrm{ref}}|$ vs. `h`. The energy involves gradient terms only (no Laplacian / biharmonic), so the inherited convergence order is the gradient's `O(h²)` (Phase 2 default). The roadmap bar "matches the analytical value" is discharged by a relative-error gate at `N = 64` (`< 1e-3` rel) plus the convergence sweep.

- **Demo binary ships, not test-only.** `examples/elastic_energy.cpp` mirrors Phase 12's CH demo and Phase 14's heterogeneous-diffusion demo — gated by `GRIDCALC_BUILD_EXAMPLES=ON` (already in CI). CLI options: `--n-x`, `--lambda`, `--mu`, `--alpha`, `--axis {x|y|z}`, `--out-dir`. Produces a single `elastic_energy.vtk` snapshot of the strain $\varepsilon_{xx}$ field via the Phase 12 `examples/common/vtk_writer.hpp`; prints the computed energy and the analytical reference. Three committed PNG snapshots in the User Guide chapter (one each at `α ∈ {0.005, 0.01, 0.02}` showing the cos-wave strain profile).

- **No new spectral cross-check entries.** Phase 9's FD–FFT fixture parameterizes over scalar-input FD operators; the Phase 15 `func::evaluate(Field<Vec3d>, callable)` surface is a higher-level wrapper, not a new differential operator. The existing rank-2 `gradient(Field<Vec3d>)` is already covered by Phase 13's identity-tests and the existing FD path; spectral coverage of vector gradients is out of Phase 15's scope.

- **No `Order` template parameter on the new `func::evaluate` overloads.** They forward to the existing `diff::gradient(Field<Vec3d>)` overload at its default order 2. A user who wants 4th-order can `materialize(...)` the gradient explicitly with `<Order=4>` and pre-compute the strain tensor — the ET layer participates either way. Adding an `Order` template parameter directly on `func::evaluate` is a one-line additive change later; no need to ship it speculatively.

- **No new build flag.** All new headers compile under the existing `gridcalc` interface target. `GRIDCALC_BUILD_EXAMPLES=ON` (already in CI) gates the new demo binary.

## Non-goals

- **No rank-3 / rank-4 tensor types.** No `core::Tensor3` / `core::Tensor4`, no `Eigen::TensorFixedSize` import. The elastic-energy demo doesn't need them.
- **No `gradient(Field<core::Mat3d>)` overload.** Rank-raising past rank-2 requires the rank-3 type that's deferred above.
- **No `func::evaluate(Field<core::Mat3d>, ...)` overload.** Same reason: the natural callable would consume `(Mat3d m, Tensor3 grad_m)` and `Tensor3` doesn't exist yet.
- **No literal rigid-stretch demo on `[0, L]³`.** Periodicity-only; periodic uniaxial wave is the substitute. Rigid stretch waits for Phase 19.
- **No nonlinear-elasticity, no plasticity, no hyperelasticity, no damage, no fracture.** Linear isotropic elasticity only — the smallest model that exercises the rank-2 contraction algebra.
- **No 4th-order constitutive stiffness $C_{ijkl}$.** Implicit in $\boldsymbol{\sigma} = \lambda \mathrm{tr}(\boldsymbol{\varepsilon})\mathbf{I} + 2\mu \boldsymbol{\varepsilon}$ for the isotropic case; no need to materialize.
- **No anisotropic elasticity** (cubic, hexagonal, triclinic). Out of scope.
- **No expression-template SIMD / OpenMP / cache-blocking.** Phase 21's performance pass.
- **No `expression to expression` arithmetic helpers beyond what the demo needs.** No vectorized `cross`, no `det`, no `inverse` on the ET layer. Add later when a use case appears.

## Deferred questions

- **Rank-3 / rank-4 tensor types.** Pinned to whichever later phase actually consumes them — likely Phase 21 (performance pass touches the high-rank surface) or a yet-to-be-added phase covering anisotropic / hyperelastic energies. Decision between in-house `Tensor3` / `Tensor4` aliases and `Eigen::TensorFixedSize` stays deferred along with the type.
- **`func::evaluate(Field<core::Mat3d>, ...)`.** Once rank-3 lands, this overload is additive. Same dispatch ladder, one new branch.
- **Symmetric-storage rank-2 type.** `Eigen::Matrix3d` stores all 9 components; for symmetric tensors like ε and σ, 6 unique components suffice. The 33 % storage savings is below the benchmark threshold for a v1.0 release.
- **`Order = 4` `func::evaluate(Field<Vec3d>, ...)` overload.** One-line additive change later.
- **Nonlinear-elasticity / hyperelasticity / Saint Venant-Kirchhoff / Neo-Hookean.** A future "Phase 15.x" or post-v1.0 phase. The current rank-2 algebra is a strict subset of what those need.

## Branch and version

- **Branch:** `2026-05-05-phase-15-tensor-functional`.
- **Version bump:** `0.14.5 → 0.15.0`.
- **CHANGELOG entry:** new `## 0.15.0 — Phase 15` block.
- **`\since` tag:** every new public symbol carries `\since 0.15.0`. Phase 13's existing `tensor::*` primitives are not modified; their `\since 0.13.0` tags stay.
