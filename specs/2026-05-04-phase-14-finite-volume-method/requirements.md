# Requirements: Phase 14 — Finite volume method (FVM)

## Roadmap reference

This phase did not exist in the original roadmap; it is **inserted** between the now-renamed Phase 13 (Vector and tensor fields, merged on 2026-05-04 as PR #17) and the previously-numbered Phase 14 (Functional evaluation for vector/tensor data, which becomes new **Phase 15**). The full mechanical renumbering is `Phase N → Phase N+1` for `N ∈ [14, 22]`, applied to `specs/roadmap.md`, `specs/STATUS.md`, and the post-Phase-13 doc chapters that already point forward to "Phase 14" (now "Phase 15"). Historical CHANGELOG entries (0.11.0, 0.12.0, 0.13.0) are left as-is — they are version-block-scoped historical records and the new 0.14.0 block calls out the renumber.

The new Phase 14 entry inserted into `specs/roadmap.md`:

> **Goal.** Conservative cell-flux discretization of the diffusion operator with heterogeneous diffusivity `D(x)`. Rounds out the solver story with the FVM foundation that the original roadmap had implicit at later phases.
>
> **Deliverables.**
> - `gridcalc::fvm` namespace and `include/gridcalc/fvm/CellLaplacian.hpp` — `fvm::cellLaplacian(psi, D)` returning $\nabla\cdot(D\nabla\psi)$ via cell-flux discretization with face-centered harmonic-mean D-averaging.
> - `solve::diffuse(psi, D_field, dt, n_steps, tag)` — new heterogeneous-D overload reusing the existing `solve::ExplicitEuler{}` / `solve::RK4{}` integrator tags. The constant-D Phase 5 / Phase 6 overload is unchanged.
> - Acceptance tests covering: agreement with the FD path on a uniform `D` field (round-off equality on uniform Cartesian grids); convergence at order 2 against an analytical heterogeneous-D solution; conservation `Σ fvm::cellLaplacian(ψ, D) = 0` to round-off on a periodic domain.
> - User Guide chapter 14 + Developer Note chapter 13.
>
> **Acceptance.** Heterogeneous-D diffusion of a manufactured solution recovers the analytical answer at order 2 in `h`; conservation holds to round-off; uniform-D path produces the same trajectory as Phase 5 / Phase 6's FD `solve::diffuse` to round-off.

## Decisions made for this feature

- **Hybrid incorporation**, not replacement. New `gridcalc::fvm` layer hosting the cell-flux Laplacian and its diffusion-solver dispatch. The Phase 1 `diff::laplacian` + Phase 5 / Phase 6 `solve::diffuse` paths stay byte-for-byte unchanged; users with constant `D` and the existing API see no change. The new heterogeneous-D path is the only place where FVM produces different physics from FD; on constant `D` the two paths agree to round-off (validated by the agreement test).

- **Heterogeneous diffusivity `D(x)` is the headline feature.** A `Field<double>` D-coefficient is the new public surface element. Without it, FVM is structurally redundant with FD on uniform Cartesian grids; with it, the new path delivers a real new capability the existing FD operators cannot.

- **Face-centered D averaged via harmonic mean.** For each face between cells `i` and `i+1`,
  $$D_{\text{face}} = \frac{2\,D_i\,D_{i+1}}{D_i + D_{i+1}}.$$
  Standard FVM choice for diffusion (LeVeque 2007 §4.2; Patankar 1980 §4.2). Recovers the correct steady-state flux when `D` jumps discontinuously between cells, unlike the arithmetic-mean alternative which biases toward the higher-D side. Harmonic-mean is undefined when either `D_i` or `D_{i+1}` is zero; the implementation asserts `D > 0` everywhere at the cell-Laplacian entry point and surfaces this in the Doxygen contract.

- **Cell-centered values, periodic-only.** Sample placement matches the existing Phase 1 `Field<double>` storage (`(i, j, k)` sits at Cartesian `(i·hx, j·hy, k·hz)`). FVM is a discretization choice, not a re-sampling choice; the underlying `Grid` and `Field` types are unchanged. Boundary conditions are periodic only — Phase 19 (the renumbered non-periodic-BC phase) handles Dirichlet / Neumann walls and now leverages the FVM foundation as the natural starting point.

- **Diffusion solver: new heterogeneous-D overload, same tags.** `solve::diffuse(psi, D_field, dt, n_steps, tag)` with `tag ∈ {ExplicitEuler{}, RK4{}}` reusing the Phase 5 / Phase 6 integrator-tag set. Internally the heterogeneous-D path calls `fvm::cellLaplacian(psi, D)` for the RHS evaluation; the integrator tag dispatches to the existing `solve::integrate(...)` machinery. Zero new tags, zero new integrator classes.

- **CFL bound generalized.** The von-Neumann CFL check on the heterogeneous-D path replaces `D · dt · Σ_a (1/h_a²) ≤ tag.diffusionCFLLimit` with `D_max · dt · Σ_a (1/h_a²) ≤ tag.diffusionCFLLimit`, where `D_max` is the spatial max of `D`. Stricter than the constant-D bound when `D` has high spots; equal to the constant-D bound when `D` is uniform. Documented on the new overload's Doxygen block.

- **`fvm::cellLaplacian` returns the time derivative `∂_t ψ`, not just `D·∇²ψ`.** Concretely, it returns `∇·(D ∇ψ)` so that the diffusion solver's RHS callable is just `[&](const auto& y){ return fvm::cellLaplacian(y, D); }`. No factor-of-D ambiguity at the integrator boundary.

- **`Order` template parameter intentionally absent.** The Phase 7 4th-order stencil family is FD-specific (asymmetric coefficient widths around cell centers); the FVM cell-flux discretization at Phase 14 ships only the 2nd-order central-difference face-flux variant. Higher-order FVM (e.g., MUSCL reconstructions) is a separate body of work and would not be motivated by the current roadmap. The `cellLaplacian` signature is therefore a plain function, not a function template; this also keeps the heterogeneous-D path simpler to read.

- **No new spectral cross-check entries.** Phase 9's FD–FFT fixture parameterizes over scalar-input FD operators; `fvm::cellLaplacian` is a different beast (heterogeneous D, conservative form). Verifying it against a spectral reference would require a manufactured solution where both `D(x)` and `ψ(x)` are band-limited and the analytical `∇·(D ∇ψ)` is computable in closed form — doable, but the manufactured-solution acceptance test in this phase already covers that case.

- **No new build-flag.** `gridcalc::fvm` is a header-only addition behind no conditional. Always built, always tested.

## Non-goals

- **No FD replacement.** Phase 1's `diff::laplacian` and Phase 5 / Phase 6's `solve::diffuse(constant-D)` are unchanged.
- **No higher-order FVM.** 2nd-order central-difference face fluxes only.
- **No rank-2 / vector FVM.** FVM for vector-valued diffusion (e.g., heat equation on a `Field<Vec3d>`) is not in scope. Phase 13's vector / tensor surface stays separate.
- **No non-periodic BCs.** Phase 19 (the renumbered non-periodic-BC phase) covers Dirichlet / Neumann / no-flux walls and consumes the FVM foundation.
- **No implicit / IMEX time stepping.** Phase 20 (renumbered) covers implicit schemes.
- **No new public type for `D`.** `D` is a plain `Field<double>` matching the existing `psi` shape; no `DiffusionCoefficient` or `LinearOperator` wrapper.
- **No CH demo update.** Phase 12's CH demo stays scalar with constant `M`; a heterogeneous-mobility CH variant could be a separate later phase but is not motivated by Phase 14.

## Deferred questions

- **Vector / tensor FVM.** Possibly Phase 21+ if a use case appears.
- **Higher-order face-flux reconstructions** (MUSCL, WENO). No use case currently. Out of the v1.0 roadmap.
- **Heterogeneous mobility for CH.** Could land as an `--D-field` extension to `examples/cahn_hilliard.cpp`. Tracked as a possible follow-up; not committed.
- **FVM-flavored `divergence` and `gradient` operators.** Those are differential-operator surfaces with no natural FVM analog (`divergence` is rank-lowering; FVM's flux-divergence formulation works on rank-0 inputs only). Stay FD-only.

## Renumbering bookkeeping

Phase numbers ≥ 14 in the original roadmap shift up by one. Mechanical renumber applied in Group 2 to:

- `specs/roadmap.md` — primary phase list.
- `specs/STATUS.md` — "Next action" wording, "Decisions worth knowing" cross-references, "Phase progress" table.
- `docs/user-guide/chapters/13-vector-tensor-fields.tex` — "Phase 14" references → "Phase 15".
- `docs/developer-note/chapters/12-vector-tensor-fields.tex` — "Phase 14" references → "Phase 15".
- `test/integrate_test.cpp` — "Phase 20" comment → "Phase 21".
- `docs/user-guide/chapters/15-diamond-lattice.tex` — "Phase 16" placeholder reference → "Phase 17" if present.

Historical CHANGELOG entries (0.11.0, 0.12.0, 0.13.0) and prior spec dirs (`specs/2026-05-*-phase-1[123]-*/`) are **not** renumbered — they are git-immutable records of decisions made at the time. The new 0.14.0 CHANGELOG block calls out the renumber explicitly.
