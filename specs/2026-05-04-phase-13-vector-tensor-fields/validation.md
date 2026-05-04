# Validation: Phase 13 — Vector and tensor fields

## Build

- [ ] CI green on Ubuntu GCC and Ubuntu Clang (Apple Clang + MSVC remain Phase 21).
- [ ] CI's `docs-build` and `docs-lint` jobs are green; new chapters render with all internal cross-references resolving.
- [ ] `\since`-tag lint passes — every new public symbol (`core::Mat3d`, the gradient/divergence overloads, the three contraction functions) carries a `\since 0.13.0` tag.
- [ ] Phase 9's FD–FFT cross-check fixture still passes — Phase 13's new operators are vector-input variants of existing scalar operators; the fixture parameterizes over scalar inputs so no extension is required.
- [ ] Warnings clean under `-Wall -Wextra -Wpedantic -Wconversion`.
- [ ] `clang-format` and `clang-tidy` pass.

## Tests

- [ ] `Field<Mat3d>` instantiates, default-constructs to zero matrices, round-trips through the i-fastest layout (`field_test.cpp` extension).
- [ ] `MatGradientOnLinearField` — `v = (a x, b y, c z)` ⇒ rank-2 gradient is `diag(a, b, c)` to round-off.
- [ ] `MatGradientIndexConvention` — `M(0, 1) = ∂_y v_x` for `v = (sin(y), 0, 0)`.
- [ ] `MatGradientConvergence_Order2` — `O(h²)` convergence slope on `N ∈ {16, 32, 64}`.
- [ ] `MatGradientOrder4` — `Order=4` overload selects the Phase 7 stencil and shows the expected accuracy improvement at fixed `N`.
- [ ] `MatDivergenceOnLinearField` — closed-form check at three interior points.
- [ ] `MatDivergenceConvergence_Order2` — `O(h²)` convergence slope.
- [ ] `MatDivergenceOrder4` — sanity check at `Order=4`.
- [ ] `TraceOfIdentity` — `trace(I) = 3` everywhere.
- [ ] `SingleContractWithIdentity` — `A · I = A` pointwise.
- [ ] `DoubleContractWithSelf` — `A : A = ‖A‖_F²` pointwise.
- [ ] `DoubleContractIsTraceOfTransposeProduct` — `A : B = trace(Aᵀ B)` pointwise.
- [ ] `HessianIsSymmetric` — **roadmap acceptance** for the `∇×∇φ = 0` family. `gradient(gradient(phi))` matrix is symmetric within tolerance on a trig manufactured solution.
- [ ] `TraceOfGradientEqualsDivergence` — second roadmap-acceptance identity. `trace(gradient(v)) = divergence(v)` to float-equality tolerance.
- [ ] `DivergenceOfGradientEqualsLaplacian` — third roadmap-acceptance identity. `divergence(gradient(v))` matches the analytical vector Laplacian (componentwise) at order 2.
- [ ] All previously-passing tests remain green: 223 prior on `clang-debug`, 147 prior on `clang-debug-nofft`. Phase 13 adds ~14 new tests.

## Documentation

- [ ] User Guide chapter 13 (`docs/user-guide/chapters/13-vector-tensor-fields.tex`) — full chapter with motivation, new types, gradient/divergence overloads with index convention, contraction primitives, three identities with code snippets, "what this does not do" section. Replaces the legacy `13-diamond-lattice.tex` placeholder (moved to `15-diamond-lattice.tex`).
- [ ] Developer Note chapter 12 (`docs/developer-note/chapters/12-vector-tensor-fields.tex`) — five-section structure (Theory · Math derivation · Algorithm · Design decisions · References). References section non-empty (≥ 1 external citation: standard continuum-mechanics text such as Gurtin / Truesdell / Holzapfel; Eigen's `Matrix3d` reference; ISO C++ language standard for the function-overloading semantics).
- [ ] `docs/user-guide/main.tex` and `docs/developer-note/main.tex` `\input{}` lines updated; legacy `13-diamond-lattice.tex` reference replaced.
- [ ] Both LaTeX skeletons rebuild cleanly via `scripts/build-docs.sh`.
- [ ] `CHANGELOG.md` has a new `## 0.13.0 — Phase 13` block following the established format.

## Performance

- [ ] No benchmark commitments at this phase. Phase 20 will profile the new operators alongside the rest of the diff stack. The contraction implementations are deliberately materialized; Phase 14 will revisit ETs.

## Bookkeeping

- [ ] `project(... VERSION 0.13.0)` in root `CMakeLists.txt`.
- [ ] `test/version_test.cpp` asserts `0.13.0`.
- [ ] `specs/STATUS.md` updated: Phase 13 row → Done, "Last updated" date bumped, "Next action" rewritten for Phase 14, new decisions added to "Decisions worth knowing", PR list extended.
- [ ] PR opened, every CI check green, no acceptance-bar relaxation.
