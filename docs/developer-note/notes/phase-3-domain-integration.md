# Developer Note — Phase 3: Domain integration (∫ over the grid)

> Phase metadata. Spec dir: `specs/2026-05-02-phase-3-domain-integration/`.
> Implements `gridcalc::func::integrate` with two reducer tags. Single
> commit version target: `0.3.0`.

## 1. Theory

Let $\Omega = [0, L_x) \times [0, L_y) \times [0, L_z)$ be a periodic
Cartesian domain with $L_a = N_a h_a$ for $a \in \{x, y, z\}$, and let
$f \in C^\infty_{\mathrm{per}}(\Omega)$ be a smooth periodic integrand.
Phase 1's `Grid` samples $f$ at the lattice points
$\mathbf{x}_{ijk} = (i h_x,\, j h_y,\, k h_z)$ for $0 \le i < N_x$,
$0 \le j < N_y$, $0 \le k < N_z$.

The continuous integral $\int_\Omega f\, dV$ is approximated by a
**Cartesian-product midpoint rule** in three dimensions. On equispaced
samples of a periodic function the midpoint rule and the trapezoidal rule
**coincide** — the trapezoidal rule's two endpoint half-weights add up to
one full-weight sample because periodicity identifies the right endpoint
with the left. This is why "periodic midpoint" and "periodic trapezoidal"
are interchangeable names for the same scheme on this grid.

The quadrature inherits the **spectrally accurate** behavior of the 1D
periodic trapezoidal rule: for analytic periodic integrands the error
decays exponentially in the grid spacing rather than algebraically. For
finitely smooth periodic integrands ($C^k$) the rate becomes $O(h^{k+1})$,
faster than the algebraic $O(h^2)$ that the rule achieves on a
non-periodic interval. See reference [1] for the full argument and the
equivalence to truncating Fourier series.

## 2. Math derivation

### Discrete formula

Apply the 1D periodic trapezoidal rule along each axis and Fubini:

$$
\int_\Omega f\, dV \;\approx\; h_x h_y h_z \sum_{i=0}^{N_x-1}
\sum_{j=0}^{N_y-1} \sum_{k=0}^{N_z-1} f_{ijk}
\;=\; \Delta V \cdot \sum_{n=0}^{N-1} f_n \;=:\; I_h,
$$

where $\Delta V = h_x h_y h_z$ is the cell volume (provided by
`Grid::getCellVolume()`) and $f_n$ are the same values laid out in the
i-fastest storage order from Phase 1. The summation is over $N = N_x N_y N_z$
elements.

### Convergence

For $f \in C^\infty_{\mathrm{per}}(\Omega)$ analytic in a strip
$|\mathrm{Im}\,z_a| < a$ for each axis, the 1D periodic trapezoidal error
satisfies $|I - I_h| \le C \exp(-2\pi a / h)$ along that axis (reference [1],
Theorem 3.2). By tensor product, the 3D error is bounded by the worst-axis
exponential. In practice this means **single-axis trig integrands and
products of single-axis trig integrands evaluate to their analytical values
within rounding** — the key observation behind Phase 3's
`SinOverPeriodIsZero` and `ManufacturedSolutionPeriodicProduct` tests, which
target an exact-zero analytical answer rather than a convergence sweep.

### Rounding-error bounds for the reducers

For a sum $S = \sum_{n=0}^{N-1} f_n$ of $N$ floating-point values
$f_n \in \mathbb{R}$, the IEEE-754 forward error of the **naive
left-to-right accumulator** is bounded by $|\hat S - S| \le \gamma_N |S^*|$
where $\gamma_N = N\varepsilon / (1 - N\varepsilon)$ and $|S^*| =
\sum_n |f_n|$ (reference [2], §4.2). The error therefore grows
**linearly** in $N$.

**Pairwise summation** (recursive halve-and-sum) replaces $\gamma_N$ with
$\gamma_{\log_2 N}$ (reference [2], §4.3, Theorem 4.4). For our typical
problem sizes $N = 32^3 \approx 3 \times 10^4$ to $128^3 \approx 2 \times 10^6$,
this is a $\sim 10^4 \times$ improvement in the worst-case bound at zero
extra arithmetic compared with the naive accumulator.

**Compensated summation** (Kahan, with Neumaier's improvement) replaces
$\gamma_N$ with $2\varepsilon$ — the bound is **independent of $N$** to
first order (reference [2], §4.4; the Neumaier variant in reference [3]
improves the unconditional inequality on the compensation step). The cost
is one extra subtraction, one comparison, and one addition per element —
empirically $\sim 3$–$4 \times$ the naive accumulator.

### Why both are exposed

Pairwise is the cheap default and is sufficient for typical grid-calculus
integrands, where values are smooth and bounded (within a few decades of
each other). Kahan is exposed because Phase 11+ functionals will integrate
$|\nabla \psi|^4$-style integrands, where the per-element values can span
many orders of magnitude in regions of sharp gradient — there, the
$O(\varepsilon \log n)$ pairwise bound can erode meaningful precision in
the smallest-magnitude regions. Closing the API at Phase 3 with both
reducers avoids a public-API revision later.

## 3. Algorithm

The reducer implementations live in
[`include/gridcalc/func/Integrate.hpp`](../../../include/gridcalc/func/Integrate.hpp),
under `gridcalc::func::detail`.

**`pairwiseSum(data, n)`** is a recursive split-and-sum with a
small-block base case at `kBaseCase = 64`. The base case runs a naive
accumulator over its block; for `n <= kBaseCase`, the routine reduces to
that accumulator, which is acceptable because the bound at `n = 64` is
$\gamma_{64} \sim 7\varepsilon$. The base-case threshold is large enough
to amortize the recursion's call-stack overhead and small enough that the
in-block error remains bounded; Phase 20 (performance) may revisit if
benchmarks indicate a different breakpoint dominates after vectorization.
Total work is $O(N)$; recursion depth is $O(\log N)$ and easily fits the
call stack at $N \le 256^3$.

**`neumaierSum(data, n)`** is a single-pass loop with a running
compensation `c`. At each step, both `|sum| >= |x|` and `|sum| < |x|` cases
are handled — the latter is Neumaier's improvement over the original Kahan
formulation (which assumed `|sum| >= |x|` and silently lost precision on
the reverse). The compensation is added back into `sum` once at the end.

**`integrate()`** is two `inline` overloads, dispatched on the empty tag
struct (`Pairwise{}` is the default-constructed default; `Kahan{}` is the
explicit second tag). Each overload reads `Grid::getCellVolume()` once and
multiplies the reducer's output by it.

**Tests** in
[`test/integrate_test.cpp`](../../../test/integrate_test.cpp). The
`DeterministicAcrossInvocations` test asserts bit-equality across two
back-to-back invocations of each reducer. At Phase 3 this is trivially
satisfied — single-threaded execution produces deterministic output by
definition. The test exists as a placeholder for Phase 20's parallel
implementation, where the same equality must hold while
`OMP_NUM_THREADS` varies; replacing the body of this test is the primary
Phase 20 deliverable for `func::integrate`.

## 4. Design decisions

- **Both `Pairwise` and `Kahan` exposed.** Pairwise as default for cheap
  use; Kahan reserved for high-dynamic-range integrands anticipated in
  Phase 11+ functionals. See the
  [`requirements.md`](../../../specs/2026-05-02-phase-3-domain-integration/requirements.md)
  "Decisions made for this feature" entry on reduction strategy.
- **Single-threaded at Phase 3.** Pulls Phase 20 work earlier if not.
  The cross-thread-count determinism contract is satisfied trivially with
  one thread; the placeholder test documents the intent. See
  `requirements.md`'s threading decision.
- **Field-only API.** `integrate(Field<Vec3d>)` and the callable
  `integrate(Field, F)` overload are deferred to Phase 11+ and Phase 4
  respectively. Anticipating them now would add API surface that may
  require revision once those phases pin down their actual signatures.

## 5. References

[1] Trefethen, L. N., & Weideman, J. A. C. (2014). "The Exponentially
Convergent Trapezoidal Rule." *SIAM Review*, 56(3), 385-458.
<https://doi.org/10.1137/130932132>. Theorem 3.2 (exponential convergence
for analytic periodic integrands); Section 4 (extension to higher
dimensions via tensor product).

[2] Higham, N. J. (2002). *Accuracy and Stability of Numerical Algorithms*
(2nd ed.). SIAM. ISBN 978-0-89871-521-7.
<https://doi.org/10.1137/1.9780898718027>. §4.2 (naive summation forward
error), §4.3 (pairwise summation, Theorem 4.4), §4.4 (compensated
summation, Kahan and improvements).

[3] Neumaier, A. (1974). "Rundungsfehleranalyse einiger Verfahren zur
Summation endlicher Summen." *Zeitschrift für angewandte Mathematik und
Mechanik*, 54(1), 39-51.
<https://doi.org/10.1002/zamm.19740540106>. Original derivation of the
improved compensation step used in `detail::neumaierSum`.

[4] Wikipedia: *Kahan summation algorithm*.
<https://en.wikipedia.org/wiki/Kahan_summation_algorithm> (permanent URL).
Accessible reference summarizing both Kahan and Neumaier variants in one
place.
