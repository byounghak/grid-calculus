# Requirements: Fix silent Grid-mismatch on RHS return in time integrators

## Bug report

`solve::integrate(...)` (both `RK4` and `ExplicitEuler` overloads) calls
the user-supplied `rhs(psi)` callable each stage and then loops
`n = 0..N-1` (with `N = psi.getSize()`) over the returned field's
`data()` pointer:

- `solve/RK4.hpp:116-125` builds `k_i = rhs(...)` four times per step;
  `solve/RK4.hpp:128-135` does the final accumulation.
- `solve/RK4.hpp:43-55` (`detail::axpyFresh`) loops `pa[n] + scale * pb[n]`.
- `solve/ExplicitEuler.hpp:97-105` does `tmp = rhs(psi); psi[n] += dt * tmp[n]`.

The static assertion
`std::is_invocable_r_v<Field<double>, Rhs&, const Field<double>&>`
only constrains the *return type*, not the runtime `Grid` carried by
the returned field. If the RHS callable returns a `Field<double>`
allocated against a different `Grid` than `psi`, four failure modes:

1. **OOB read** if any returned `k_i` is allocated on a grid with
   *fewer* total cells than `psi.getSize()`. The `pb[n]` /
   `tmp.data()[n]` reads pass the end of the smaller field's
   `std::vector` storage. Undefined behavior.
2. **Silent layout drift** if the returned field has the same total
   size but different per-axis extents (e.g., `psi` is `8×8×8`, `k_i`
   is `4×16×8`; both have 512 elements). No OOB, but the i-fastest
   linear index `i + Nx * (j + Ny * k)` interprets the same flat `n`
   as physically different sample points across `psi` and `k_i`. The
   accumulation silently mixes derivatives from misaligned cells.
3. **Silent `1/h` drift** if `psi.getGrid()` and the RHS's captured
   grid have the same shape but different cell sizes. The RHS's
   spatial derivatives were scaled by the wrong `1/h`; sizes match,
   no OOB, no warning — `psi` accumulates wrong-magnitude RHS values.
4. **Silent reduction with another field on a different
   discretization** if the RHS callable couples to a captured second
   field (e.g., heterogeneous-D `D_field`) sized to a different grid.
   `axpyFresh` happily folds them together at the same flat `n`.

`detail::axpyFresh` is the lower-level chokepoint shared by all three
RK4 stage builders; its docstring says "must have the same shape" but
enforces nothing.

## Why the codebase is unaffected today (dormant bug)

Every RHS that the project's own operators produce is grid-preserving
by construction:

- `diff::laplacian`, `gradient`, `divergence`, `biharmonic`, every
  Phase 8 d-prefix partial — all return `Field<U>(field.getGrid())`
  and propagate the grid.
- `fvm::cellLaplacian(psi, D)` returns a fresh field on
  `psi.getGrid()`.
- The repo's RHS usages (`test/cahn_hilliard_test.cpp`,
  `test/heterogeneous_diffusion_test.cpp`,
  `test/diffusion_test.cpp`, `test/integrator_test.cpp`,
  `examples/heterogeneous_diffusion.cpp`,
  `examples/cahn_hilliard.cpp`) all build their RHS as a closure over
  `psi`'s grid via the diff:: / fvm:: operators, so the returned grid
  matches by construction.

To trigger the bug, user code would have to write an RHS that:

- Holds a captured `Grid` different from `psi.getGrid()`; or
- Allocates against a singleton / global `Grid`; or
- Couples to another field on a different discretization.

The project's mission target ("production / industrial") is exactly
the audience that *will* eventually write such code. Catching this at
the integrator boundary is cheap; debugging a corrupted state vector
mid-run is not.

## Decisions made for this feature

- **Failure mode: `throw std::invalid_argument`** at the integrator
  entry per stage and at `detail::axpyFresh` entry, when the consumed
  field's `Grid` does not bit-match. Matches the existing `Grid`-
  constructor pattern (`include/gridcalc/core/Grid.hpp`); fires
  identically in Debug and Release; the integrator is a hot loop so
  the check is one branch per stage, not per grid point.

- **Check location: integrator entry + `axpyFresh` (defense-in-depth).**
  - `axpyFresh` is the genuine low-level chokepoint — every grid-
    mismatched fold-in must go through it (or the integrator's final
    accumulation, which gets its own check).
  - Each integrator's per-stage `rhs(...)` call gets a check
    immediately after, with an error message naming the offending
    stage (`k1` / `k2` / `k3` / `k4` / `tmp`) so the user can locate
    the bug from the message alone.
  - The final `psi += sum(k_i)` accumulation in `RK4::integrate` gets
    one check at the start (before any stage runs); axpyFresh covers
    the per-stage builds.

- **`Grid` equality is bit-exact, componentwise.** New
  `Grid::operator==(const Grid&) const noexcept` and `operator!=`.
  Compares `Nx`, `Ny`, `Nz`, and the three components of `cell_size`
  with `==`. Rationale:
  - An RHS produced by composing `diff::*` / `fvm::*` operators on
    `psi.getGrid()` returns a `Field` whose `Grid` is bit-copied from
    `psi`'s; bit-exact equality holds by construction.
  - There is no source of numerical drift on `Grid` values in this
    codebase — they are constructed once at the top of a simulation
    and passed around by const reference. A non-bit-exact mismatch is
    *always* a real bug (different `Grid` object), not a numerical
    artifact, so we want to catch it. A tolerance-based comparison
    would mask exactly the class of bugs the check exists to find.

- **Narrow scope: integrator + `axpyFresh` only.** The
  diff:: operators that produce RHS values already do
  `Field<U> result(field.getGrid())` and propagate the grid by
  construction. `tensor::*` contractions and `func::evaluate` consume
  two `Field` arguments and assume their grids match — worth a
  separate audit (Phase 21 / dedicated follow-up), but not folded
  here. The 0.14.2 patch closes the integrator-side gap; broader
  audit is its own work item.

- **Shared helper.** New free function
  `gridcalc::solve::detail::requireSameGrid(const Field<double>& a, const Field<double>& b, const char* context)`
  in `include/gridcalc/solve/detail/RhsGridCheck.hpp`. Throws
  `std::invalid_argument` with a message that names the context
  (e.g., `"solve::integrate(... RK4): k2"`), the expected and actual
  `(Nx, Ny, Nz, hx, hy, hz)` tuples, so the user can diagnose without
  reading source. Reused by `axpyFresh` and by both integrator
  overloads. Lives under `solve/detail/` to keep it out of the public
  API surface (mirrors `diff/detail/PreconditionAxisExtent.hpp` from
  0.14.1).

- **Patch release `0.14.2`.** Second patch release on the `0.14`
  line. CMake `VERSION` bumps from `0.14.1` to `0.14.2`; CHANGELOG
  gains a `0.14.2` block under "Fixed". The new `Grid::operator==`
  carries `\since 0.14.2`; the new helper carries `\since 0.14.2`;
  touched integrator Doxygen blocks gain a `\throws` line.

## Non-goals

- **No tensor / func::evaluate audit.** Out of scope for this patch;
  see "Decisions made for this feature."
- **No `Grid` storage refactor.** The `_Nx, _Ny, _Nz, _cell_size`
  members stay private; equality is computed via the existing public
  getters (`getNx`, `getNy`, `getNz`, `getCellSize`). No friend
  declaration needed.
- **No integrator persistence / scratch-pool refactor.** Phase 21 is
  the natural home for switching from per-stage allocation to a
  reusable scratch pool; this patch leaves the existing six-Field-
  per-step pattern untouched and adds the precondition checks at the
  same call sites.
- **No `detail::axpyFresh` rename.** The "Fresh" suffix flags the
  allocator behavior; it stays. Only its body gains the precondition.
- **No new public type for "grid mismatch."** `std::invalid_argument`
  with a precise message matches the existing project pattern; a
  custom `GridMismatchError` exception type would force callers to
  add a second `catch` in their integrator-driver code with no real
  benefit.
- **No spectral-path check.** `gridcalc::spectral::*` is verification-
  only and does not participate in time integration.

## Deferred questions

- **`tensor::singleContract` / `tensor::doubleContract` /
  `func::evaluate` Grid-mismatch audit.** Same pattern, different
  failure mode (no time-stepping accumulation, but still silent
  flat-index reduction across mismatched grids). Tracked as a
  follow-up; folding here would expand the patch surface beyond what
  the integrator bug motivates.
- **`Field::operator==` / structural equality on `Field<T>`.** Out of
  scope; the bug is about `Grid`, not field values.
