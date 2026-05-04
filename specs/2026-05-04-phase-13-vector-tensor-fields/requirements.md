# Requirements: Phase 13 — Vector and tensor fields

## Roadmap reference

> **Goal.** Generalize `Field<T>` beyond scalars; rank-raising gradients.
>
> **Deliverables.**
> - `Field<Vector3d>`, `Field<Matrix3d>` instantiations.
> - `gradient(Field<Vector3d>) -> Field<Matrix3d>` (∂_i v_j).
> - `divergence(Field<Matrix3d>) -> Field<Vector3d>`.
> - `tensor/Contraction.hpp` — single and double contractions.
>
> **Acceptance.** Vector identities (e.g. $\nabla\times\nabla\phi = 0$) recovered to discretization error.

(`specs/roadmap.md`, Phase 13.)

## Decisions made for this feature

- **Rank-2 type: `Eigen::Matrix3d` directly, exposed as `gridcalc::core::Mat3d` in `core/EigenAliases.hpp`.** Standard Eigen module, not the unsupported tensor module — Phase 13 stays within the dependency surface already pinned at Phase 0. The `Mat3d` alias lands alongside `Vec3d` (Phase 1 / 2) per the existing decision recorded under "Decisions worth knowing". Rank-3 / rank-4 tensor types are **not** introduced in Phase 13; they are deferred to Phase 14 alongside the contraction expression-template work that motivates them. The `TensorFixedSize` portability risk in `STATUS.md`'s "Open / deferred items" stays deferred.

- **Index convention for the rank-2 gradient: `M(i, j) = ∂_j v_i`.** Row index `i` selects the vector component; column index `j` selects the derivative axis. Trace gives the divergence (`tr(M) = ∂_x v_x + ∂_y v_y + ∂_z v_z = ∇·v`) under either ordering, but this convention is the standard CFD / continuum-mechanics velocity-gradient tensor `L = ∇v`. Documented in the Doxygen block on the new gradient overload and in the User Guide chapter.

- **Index convention for the rank-2 divergence: `(div M)_i = ∂_j M_ij`** (sum over the second index). With the rank-2 gradient convention above, this immediately gives `div(grad(v))_i = ∂_j ∂_j v_i = (vector Laplacian)_i` — the third acceptance identity drops out by construction.

- **Naming: function overloading on input type.** `diff::gradient` and `diff::divergence` get new overloads on `Field<core::Vec3d>` / `Field<core::Mat3d>` respectively. Same names, dispatched by input type. Phase 14 may add a unifying templated form when rank-3 inputs arrive; for Phase 13's two cases, plain function overloading is the smallest delta against the Phase 2 surface. The Phase 7 `Order` template parameter is preserved on both overloads (default `2`, `4` available); the overloads share the same stencil tables as the scalar versions.

- **Contractions: fully-materialized in Phase 13; expression templates revisited in Phase 14.** Mirrors the Phase 11 deferral rationale. `tensor/Contraction.hpp` ships:
  - `trace(Field<Mat3d>) -> Field<double>` — sum of diagonals; the ergonomic helper for the trace=divergence identity test.
  - `singleContract(Field<Mat3d>, Field<Mat3d>) -> Field<Mat3d>` — pointwise matrix product (`(A·B)(i,j) = A(i,k) B(k,j)`).
  - `doubleContract(Field<Mat3d>, Field<Mat3d>) -> Field<double>` — pointwise full contraction (`A:B = A(i,j) B(i,j)`); equals `tr(A·Bᵀ)`. The natural choice for elastic energy density `½σ:ε` in Phase 14.
  
  All three return fresh `Field<...>` outputs. Rank-2 inputs only; the rank-6 functional-evaluation explosion that motivates ETs does not yet appear, so the materialized form is appropriate. Phase 14's spec round will add ETs alongside `func::evaluate(Field<Vec3d>, ...)`.

- **Acceptance bar: three vector identities, programmatic in CI.** Each is a one-or-two-test convergence sweep on a manufactured field where the analytical answer is known:
  1. **Hessian symmetry** (the literal roadmap example, $\nabla\times\nabla\phi = 0$, expressed in rank-raising vocabulary). For $\phi$ smooth periodic, $H := \nabla\nabla\phi$ is symmetric: $H(i,j) = H(j,i)$. Test: `gradient(gradient(phi))` is symmetric within tolerance for a trig manufactured solution.
  2. **Trace = divergence.** $\mathrm{tr}(\nabla v) = \nabla\cdot v$. Test on a manufactured vector field where both sides are computable in closed form.
  3. **Vector Laplacian.** $\nabla\cdot(\nabla v) = \nabla^2 v$ (componentwise). Test against `diff::laplacian` applied to each Cartesian component of the same manufactured field.
  
  Plus per-operator order-2 convergence sweeps (gradient on `Vec3d`, divergence on `Mat3d`, trace, single/double contractions) to discharge the standing accuracy guarantee.

## Non-goals

- **No `Field<Tensor3>` or `Field<Tensor4>` of any flavor.** Deferred to Phase 14.
- **No expression-template contraction layer.** Deferred to Phase 14.
- **No curl operator as a free function.** Phase 13's `gradient(gradient(phi))` symmetry test discharges the $\nabla\times\nabla\phi = 0$ identity in rank-raising vocabulary; an explicit `curl` (rank-1 → rank-1) is not on the roadmap and is not added speculatively.
- **No vector-valued IndexPolicy work.** Phase 18's Dirichlet / Neumann work covers boundary conditions; for Phase 13, the periodic policy on `Field<Mat3d>` and `Field<Vec3d>` is the only one that exists.
- **No new spectral cross-check entries.** The Phase 9 fixture parameterizes over scalar-input FD operators; the rank-raising operators have no Phase 9 analog yet (spectral derivatives of vector fields would be a separate, Phase 14-or-later body of work).
- **No CH-or-other example update.** Phase 12's CH demo stays scalar; Phase 14 introduces the first vector/tensor demo (linear elastic energy).
- **No 4th-order stencil rollout for the new operators.** The new gradient and divergence overloads accept the `Order` template parameter (default 2, 4 available) but only `Order=2` is exercised in the Phase 13 acceptance tests. A 4th-order convergence sweep on each new operator could be added later if a use case appears.

## Deferred questions

- **Rank-3 / rank-4 tensor types and their contractions.** Phase 14 — the natural home alongside the func::evaluate-over-tensor-fields surface that finally needs them. The choice between in-house aliases vs. Eigen unsupported `TensorFixedSize` reopens at that point.
- **Contraction expression templates.** Phase 14 — see above.
- **Spectral cross-check for vector / tensor operators.** Possibly Phase 14 if it benefits from a manufactured-solution gate; otherwise deferred.
- **Documentation chapter renumbering.** The existing placeholder `docs/user-guide/chapters/13-diamond-lattice.tex` (legacy from Phase 11's renumber) is moved to `15-diamond-lattice.tex` to free chapter slot 13 for this phase. Phase 14 will land `14-tensor-functional.tex` between them.
