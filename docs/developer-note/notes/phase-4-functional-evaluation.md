# Developer Note — Phase 4: Functional evaluation, callable API

> Phase metadata. Spec dir:
> `specs/2026-05-02-phase-4-functional-evaluation/`. Implements
> `gridcalc::func::evaluate`. Version target: `0.4.0`.

## 1. Theory

A **functional** maps a function (here, a scalar field) to a real number:

$$
F[\psi] \;=\; \int_{\Omega} f\bigl(\psi(\mathbf{x}),\, \nabla\psi(\mathbf{x}),\, \nabla^2\psi(\mathbf{x})\bigr)\, dV,
$$

with $\psi : \Omega \to \mathbb{R}$ on the periodic Cartesian domain
$\Omega = [0, L_x) \times [0, L_y) \times [0, L_z)$ and $f$ a smooth
density. Functionals of this form arise throughout continuum and
statistical physics: Ginzburg–Landau and Cahn–Hilliard free energies, the
Helfrich bending energy of membranes, kinetic-energy functionals in
density-functional theory, and many phase-field models. Phase 4 stands up
the **scalar-field, periodic** version; vector / tensor / non-periodic
extensions arrive in later phases.

The **discrete functional** $F_h$ replaces the integral with the periodic
midpoint quadrature from Phase 3 and the continuous derivatives with the
2nd-order central differences from Phases 1 and 2:

$$
F_h(\psi) \;=\; \Delta V \sum_{i,j,k} f\bigl(\psi_{ijk},\, (\nabla_h \psi)_{ijk},\, (\nabla^2_h \psi)_{ijk}\bigr).
$$

The convergence of $F_h \to F$ as $h \to 0$ is bounded by the **slowest
component**. The integration is spectrally accurate for analytic periodic
integrands (Phase 3, reference [3]). The gradient and Laplacian stencils
are 2nd-order. The composed convergence rate is therefore $O(h^2)$
**whenever the integrand actually consumes a derivative**; integrands that
depend on $\psi$ alone inherit the spectrally accurate quadrature error
directly.

## 2. Math derivation

### Convergence

For a smooth functional density $f$, Taylor-expand around the exact
derivatives:

$$
f(\psi, \nabla_h\psi, \nabla^2_h\psi) \;=\; f(\psi, \nabla\psi, \nabla^2\psi) \;+\; \partial_{(\nabla\psi)} f \cdot (\nabla_h - \nabla)\psi \;+\; \partial_{(\nabla^2\psi)} f \cdot (\nabla^2_h - \nabla^2)\psi \;+\; \text{higher order}.
$$

The pointwise gradient error is $(\nabla_h - \nabla)\psi = O(h^2)$ (Phase
2 stencil); the pointwise Laplacian error is $O(h^2)$ (Phase 1 stencil);
and integrating an $O(h^2)$ pointwise error over a finite domain produces
an $O(h^2)$ functional error (the integration is spectrally accurate for
the smooth periodic prefactor). For the Ginzburg–Landau example below,
both the leading-order constant and the $O(h^2)$ convergence are matched
by experiment to within 5%.

### Hand-computed reference for Ginzburg–Landau

For $\psi = \sin(x)\sin(y)\sin(z)$ on $[0, 2\pi]^3$ and $W = 1$:

**Gradient term** $\int (\tfrac{1}{2})|\nabla\psi|^2\, dV$:

$$
|\nabla\psi|^2 \;=\; \cos^2 x \sin^2 y \sin^2 z \;+\; \sin^2 x \cos^2 y \sin^2 z \;+\; \sin^2 x \sin^2 y \cos^2 z.
$$

Using $\int_0^{2\pi}\sin^2 = \int_0^{2\pi}\cos^2 = \pi$, each
axis-aligned $\cos^2\sin^2\sin^2$ tensor product integrates to
$\pi \cdot \pi \cdot \pi = \pi^3$. Three terms total give
$3\pi^3$, hence

$$
\int \tfrac{1}{2}|\nabla\psi|^2\, dV \;=\; \tfrac{3}{2}\pi^3.
$$

**Potential term** $\int (\psi^2 - 1)^2\, dV = \int [\psi^4 - 2\psi^2 + 1]\, dV$:

- $\int \psi^4\, dV$ uses $\int_0^{2\pi}\sin^4 = 3\pi/4$, so
  $\int \sin^4 x \sin^4 y \sin^4 z\, dV = (3\pi/4)^3 = 27\pi^3/64$.
- $\int \psi^2\, dV = \pi \cdot \pi \cdot \pi = \pi^3$.
- $\int 1\, dV = (2\pi)^3 = 8\pi^3$.
- Sum: $27\pi^3/64 - 2\pi^3 + 8\pi^3 = 27\pi^3/64 + 6\pi^3 = 411\pi^3/64$.

**Total** (with $W = 1$):

$$
F = \tfrac{3}{2}\pi^3 \;+\; \tfrac{411}{64}\pi^3 \;=\; \tfrac{507}{64}\pi^3 \;\approx\; 245.62.
$$

### Why the 1-arg path is exact

When $f$ depends only on $\psi$, no stencil is applied — $f(\psi_{ijk})$
is read directly from the input field, and the only error comes from the
Phase 3 periodic-trapezoidal rule. For analytic periodic integrands the
quadrature converges exponentially in $h$ (reference [3], Theorem 3.2),
so the test
[`FunctionalTest.PsiOnlyArity`](../../../test/functional_test.cpp) sets
its tolerance to $10^{-12}$ rather than the $10^{-2}$ used by the
gradient-bearing tests.

## 3. Algorithm

The implementation lives in
[`include/gridcalc/func/Functional.hpp`](../../../include/gridcalc/func/Functional.hpp).

```cpp
template <class F, class Tag = Pairwise>
inline double evaluate(const core::Field<double>& psi, F&& f, Tag tag = {});
```

### Arity dispatch

The key piece is an `if constexpr` ladder on
`std::is_invocable_r_v<double, F&, …>`:

1. `f(double, const Vec3d&, double)` matches → materialize **both**
   `diff::gradient(psi)` and `diff::laplacian(psi)`.
2. else `f(double, const Vec3d&)` matches → materialize **only**
   `diff::gradient(psi)`.
3. else `f(double)` matches → no stencil materialization.
4. else trip a `static_assert` via the `detail::deferred_false<F>` trick.

Listing the longest-arity branch first means a callable with default
arguments (which would match all three signatures) selects the full
triplet, which is the safer interpretation — the user provided defaults
specifically to allow the extra arguments.

The `detail::deferred_false<F>` indirection is required because plain
`static_assert(false, …)` inside an `if constexpr` else-branch fires
unconditionally under C++17 (resolved in C++23 by CWG2518; reference
[4]).

### Eager materialization

Each arity branch:

1. Calls `diff::gradient(psi)` and/or `diff::laplacian(psi)` as needed.
   These return fresh `Field<Vec3d>` and `Field<double>` respectively, so
   total temporary storage is up to $4 N_x N_y N_z$ doubles (3 for the
   gradient field's Vec3d entries, 1 for the Laplacian).
2. Allocates a `Field<double> integrand(grid)` of the same size.
3. Iterates `(k, j, i)` (z-slowest, x-fastest — matches the i-fastest
   storage layout from Phase 1) and writes `f(...)` into `integrand`.
4. Returns `func::integrate(integrand, tag)`.

Algorithmic complexity is $O(N_x N_y N_z)$. Memory is $O(N_x N_y N_z)$
on top of the input.

### Tag forwarding

`evaluate` takes a `Tag` defaulted to `Pairwise{}` and forwards it to
`integrate`. A `static_assert` at the top of the function constrains
`Tag` to `Pairwise` or `Kahan` so a misspelled tag fails at the call
site with a readable diagnostic instead of an unrelated overload-set
error from `integrate`.

### Tests

`test/functional_test.cpp` exercises all three dispatch paths plus the
Ginzburg–Landau hand-computed reference and a 2nd-order convergence
sweep.

## 4. Design decisions

- **SFINAE-detected arity rather than the strict roadmap signature.** The
  roadmap specifies `f(double, Vector3d, double)` only; Phase 4 widens
  this to a strict superset (1-, 2-, or 3-arg) because the dispatch is
  trivial in `if constexpr`, the user gets ergonomically natural
  signatures, and unused stencils are skipped — a real correctness and
  performance win even at Phase 4. See the
  [`requirements.md`](../../../specs/2026-05-02-phase-4-functional-evaluation/requirements.md)
  "Decisions made for this feature" entry.
- **Eager materialization over per-point lazy evaluation.** Per-point
  lazy stencils would save the gradient and Laplacian field allocations
  but duplicate the stencil formulas inline, blocking direct reuse of
  `diff::gradient` / `diff::laplacian`. Phase 11+ has a stronger reason
  to revisit this (rank-extended functionals where the rank-6 Hessian of
  a rank-2 field would explode memory if materialized).
- **Vec3d argument by `const` reference.** Per-point Vec3d allocation is
  cheap on the stack but copying the 24 bytes is wasteful in the inner
  loop. The SFINAE check uses `const Vec3d&` so user lambdas can declare
  it either way (by-value lambdas accept the const-ref argument and
  copy).
- **Worked example is Ginzburg–Landau** rather than a simpler integrand:
  matches the roadmap deliverable and exercises the 2-arg dispatch path,
  which is the most common arity in practice (gradient-energy
  functionals).

## 5. References

[1] Provatas, N., & Elder, K. (2010). *Phase-Field Methods in Materials
Science and Engineering*. Wiley-VCH. ISBN 978-3-527-40747-7. Chapter 3
introduces the Ginzburg–Landau functional in the form
$F = \int (\tfrac{1}{2})|\nabla\psi|^2 + W\,(\psi^2-1)^2\, dV$ used as
the worked example here. <https://doi.org/10.1002/9783527631520>.

[2] Boyd, J. P. (2001). *Chebyshev and Fourier Spectral Methods* (2nd
ed.). Dover. ISBN 978-0-486-41183-5. Chapter 4 covers tensor-product
quadrature and the inheritance of convergence rates from the
slowest-component approximation, the argument used in §1 of this note.
Free PDF at the author's site:
<https://websites.umich.edu/~jpboyd/BOOK_Spectral2000.html> (permanent
URL).

[3] Trefethen, L. N., & Weideman, J. A. C. (2014). "The Exponentially
Convergent Trapezoidal Rule." *SIAM Review*, 56(3), 385–458.
<https://doi.org/10.1137/130932132>. Used here for the spectral
accuracy of the integration step (cited in Phase 3's Developer Note as
well).

[4] cppreference: `std::is_invocable_r`.
<https://en.cppreference.com/w/cpp/types/is_invocable> (permanent URL).
Provides the precise semantics of the SFINAE check used in the arity
dispatch.

[5] CWG Issue 2518 (C++ Core Working Group). "Conformance requirements
and `static_assert(false)` in `if constexpr`."
<https://www.open-std.org/jtc1/sc22/wg21/docs/cwg_defects.html#2518>
(permanent URL). Documents the C++23 resolution that allows
`static_assert(false, …)` inside `if constexpr` without the
`deferred_false` workaround. Phase 4 still targets C++17 per
`tech-stack.md`, so the workaround is required.
