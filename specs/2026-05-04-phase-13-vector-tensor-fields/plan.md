# Plan: Phase 13 — Vector and tensor fields

Seven commit-sized groups. Build and tests stay green between every commit.

## Group 1 — Spec files

- `specs/2026-05-04-phase-13-vector-tensor-fields/{requirements, plan, validation}.md` (this directory).

**Commit message:** `docs: phase 13 spec files (vector and tensor fields)`

## Group 2 — Mat3d alias and Field<Mat3d> instantiation

- `include/gridcalc/core/EigenAliases.hpp` — add `using Mat3d = Eigen::Matrix3d;` alongside the existing `Vec3d`. Doxygen block with `\since 0.13.0`.
- `test/field_test.cpp` — extend with `FieldMat3dTest` cases verifying `Field<Mat3d>` instantiates, default-constructs to zero matrices (Eigen default), accepts `operator()(i, j, k)` access for read/write, and round-trips through the i-fastest layout. No new file; just additional tests in the existing `field_test.cpp`.

**Commit message:** `feat(core): add Mat3d alias and Field<Mat3d> instantiation tests`

## Group 3 — Rank-raising gradient: Field<Vec3d> → Field<Mat3d>

- `include/gridcalc/diff/Gradient.hpp` — add an overload (or sibling template) `template <int Order = 2> Field<Mat3d> gradient(const Field<Vec3d>&)` that produces `M(i, j) = ∂_j v_i`. Implementation: per-component `gradient(v_i)` (Phase 2 scalar gradient applied to each Cartesian projection of `v`), arranged into the 3×3 matrix at every grid point. Doxygen block with `\since 0.13.0` documents the index convention.
- `test/vector_tensor_test.cpp` — new file, four tests for the rank-2 gradient:
  - `MatGradientOnLinearField` — `v = (a x, b y, c z)` ⇒ `M = diag(a, b, c)` exactly (within round-off, since the linear field has exact 2nd-order central-difference derivatives).
  - `MatGradientIndexConvention` — explicit check that `M(0, 1) = ∂_y v_x` for `v = (sin(y), 0, 0)`.
  - `MatGradientConvergence_Order2` — convergence sweep on a trig manufactured solution; assert `O(h²)` slope on `N ∈ {16, 32, 64}`.
  - `MatGradientOrder4` — quick check that `Order=4` overload selects the Phase 7 stencil and gives a tighter error at fixed `N`.

**Commit message:** `feat(diff): rank-raising gradient(Field<Vec3d>) -> Field<Mat3d>`

## Group 4 — Rank-lowering divergence: Field<Mat3d> → Field<Vec3d>

- `include/gridcalc/diff/Divergence.hpp` — add an overload (or sibling template) `template <int Order = 2> Field<Vec3d> divergence(const Field<Mat3d>&)` returning `(div M)_i = ∂_j M_ij`. Implementation: per-component `divergence(M_i_row)` (Phase 2 scalar divergence applied to each row of `M`), assembled into the `Vec3d` output. Doxygen block with `\since 0.13.0`.
- `test/vector_tensor_test.cpp` — three more tests:
  - `MatDivergenceOnLinearField` — closed-form check at three points.
  - `MatDivergenceConvergence_Order2` — convergence sweep on trig.
  - `MatDivergenceOrder4` — `Order=4` overload sanity.

**Commit message:** `feat(diff): rank-lowering divergence(Field<Mat3d>) -> Field<Vec3d>`

## Group 5 — Tensor contraction primitives

- New file `include/gridcalc/tensor/Contraction.hpp` — three free functions, all in `gridcalc::tensor` namespace:
  - `trace(const Field<Mat3d>&) -> Field<double>` — point-wise `M.trace()`.
  - `singleContract(const Field<Mat3d>&, const Field<Mat3d>&) -> Field<Mat3d>` — point-wise matrix product.
  - `doubleContract(const Field<Mat3d>&, const Field<Mat3d>&) -> Field<double>` — point-wise `(A.array() * B.array()).sum()`.
  
  Each carries a Doxygen block with `\since 0.13.0` and the explicit index formula.
- `test/vector_tensor_test.cpp` — four more tests:
  - `TraceOfIdentity` — `trace(I)` field equals 3 everywhere.
  - `SingleContractWithIdentity` — `A · I = A` pointwise.
  - `DoubleContractWithSelf` — `A : A = ‖A‖_F²` pointwise (Frobenius squared).
  - `DoubleContractIsTraceOfTransposeProduct` — `A : B = trace(Aᵀ B)` pointwise (consistency between the two contraction shapes).

**Commit message:** `feat(tensor): contraction primitives (trace, single, double)`

## Group 6 — Vector identities (the roadmap acceptance bar)

- `test/vector_tensor_test.cpp` — the three identity tests:
  - `HessianIsSymmetric` — `H := gradient(gradient(phi))` is symmetric within stencil-order tolerance for `phi(x,y,z) = sin(x) cos(y) sin(z)` on a 32³ grid (this is the rank-raising restatement of the roadmap's `∇×∇φ = 0` example).
  - `TraceOfGradientEqualsDivergence` — `trace(gradient(v)) = divergence(v)` for a trig vector field, within float-equality tolerance (both sides use the same scalar diff::gradient stencil so they should agree at machine precision modulo associativity).
  - `DivergenceOfGradientEqualsLaplacian` — `divergence(gradient(v)) = laplacian(v)` (componentwise, using Phase 1's `diff::laplacian` on each `v_i`); within stencil-order tolerance — the discrete `divergence(gradient(·))` round-trip differs from a direct Laplacian by the documented stride-2 phenomenon (Phase 2 STATUS line), so the test compares against the analytical vector Laplacian, not against `diff::laplacian` of the same field. (Mirrors the Phase 2 round-trip test against the analytical Laplacian, generalized to vector input.)

**Commit message:** `test(diff,tensor): three vector identities (Hessian symmetry, trace=div, vector Laplacian)`

## Group 7 — Doc deliverables and bookkeeping

- `docs/user-guide/chapters/13-vector-tensor-fields.tex` — full chapter (replacing the legacy `13-diamond-lattice.tex` placeholder which moves to `15-...`). Sections: motivation; new types; gradient and divergence overloads with the index convention; contraction primitives; the three identities, with code snippets; what this does not do.
- `docs/user-guide/chapters/15-diamond-lattice.tex` — placeholder content moved from `13-diamond-lattice.tex` (renamed; Phase 15's eventual real chapter will overwrite).
- `docs/user-guide/main.tex` — drop the old `\input{chapters/13-diamond-lattice.tex}`, add `\input{chapters/13-vector-tensor-fields.tex}` and `\input{chapters/15-diamond-lattice.tex}`.
- `docs/developer-note/chapters/12-vector-tensor-fields.tex` — new chapter, five-section structure (Theory, Math derivation, Algorithm, Design decisions, References). At least one external citation.
- `docs/developer-note/main.tex` — add `\input{chapters/12-vector-tensor-fields.tex}` after chapter 11.
- Bump `project(... VERSION 0.13.0)` in root `CMakeLists.txt`.
- Update `test/version_test.cpp` to assert `0.13.0`.
- New `## 0.13.0 — Phase 13` block in `CHANGELOG.md`.
- Refresh `specs/STATUS.md`: Phase 13 row → Done, "Last updated" date bump, "Next action" rewritten for Phase 14, new decisions added.

**Commit message:** `chore: bump version to 0.13.0 and refresh STATUS for Phase 13`

---

After Group 7: `git push -u origin 2026-05-04-phase-13-vector-tensor-fields`, open PR, watch CI per `specs/CLAUDE.md` step 6.
