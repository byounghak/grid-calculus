# Requirements: Phase 9 — Spectral verification harness

## Roadmap reference

> ## Phase 9 — Spectral verification harness
>
> - **Goal.** A second, independent path for computing derivatives on **periodic** grids — used as the gold-standard reference against which every finite-difference operator is automatically validated. The FFT path is **not** exposed as a primary differentiation engine; the FD operators remain the production code.
> - **Scope at this phase.** Orthogonal Cartesian periodic grids only. Non-orthogonal-lattice and multi-atom-basis spectral counterparts are added as deliverables of Phases 15 and 16 alongside the corresponding grid features.
> - **Deliverables.**
>   - `spectral/Fft.hpp` — thin wrapper over PocketFFT for forward/backward 3D real-to-complex FFTs on `Field<double>`.
>   - `spectral/Derivatives.hpp` — spectral counterparts of the periodic FD operators delivered through Phase 8:
>     - `spectral::gradient(field) -> Field<Vector3d>`
>     - `spectral::laplacian(field) -> Field<double>`
>     - `spectral::biharmonic(field) -> Field<double>`
>     - `spectral::mixed_partial(field, i, j) -> Field<double>`
>
>     Implemented as multiplication by $i\mathbf{k}$, $-|\mathbf{k}|^2$, $|\mathbf{k}|^4$, etc. in Fourier space, then inverse transform.
>   - **Wavenumber and aliasing discipline.** Standard convention $k_x = 2\pi n/(N_x h_x)$ for $n \in [-N/2, N/2)$. The Nyquist mode is zeroed for odd-order derivatives to avoid aliasing artifacts. Documented in the operator headers.
>   - **CI cross-check fixture.** A parameterized GoogleTest fixture that, for every periodic FD operator at every advertised accuracy order, applies the FFT version to the same input and verifies the FD–FFT discrepancy scales as $\mathcal{O}(h^p)$ where $p$ is the FD order. Becomes a permanent CI gate; future FD changes that silently break accuracy are caught here.
>   - **Build flag.** `-DGRIDCALC_USE_FFT=ON` (default). Disabling skips the spectral tests but does not affect FD code.

## Decisions made for this feature

- **Vendored PocketFFT.** `pocketfft_hdronly.h` is committed into `extern/pocketfft/` (header-only single file). No network at build time, fully self-contained. Trade-off vs `FetchContent` (which Eigen and GoogleTest use): manual bump to upgrade. Picked because PocketFFT is small, stable, and the project's mission target is "production / industrial," for which build determinism without network access is valuable. Upstream source: <https://gitlab.mpcdf.mpg.de/mtr/pocketfft>. Vendored revision is recorded in `extern/pocketfft/VERSION` (revision SHA + date) and in `tech-stack.md`.
- **Hybrid spectral API.** Implementation core is the generic `spectral::partial<int Nx, int Ny, int Nz>(const Field<double>& field) -> Field<double>` template — multi-index spectral partial, multiplies the FFT spectrum by $(ik_x)^{N_x}(ik_y)^{N_y}(ik_z)^{N_z}$ and inverse-transforms. Named convenience wrappers cover the contracted / rank-1 operators that already have non-d-prefix names in Phase 1, 2, and 8: `spectral::laplacian(field)`, `spectral::biharmonic(field)`, `spectral::gradient(field) -> Field<Vec3d>`. The roadmap's literal `spectral::mixed_partial(field, i, j)` is dropped in favor of the d-prefix-compatible `spectral::partial<Nx, Ny, Nz>(field)` — for $\partial^2/\partial x \partial y$ that is `spectral::partial<1, 1, 0>(field)`. **No 31 named `spectral::*` d-prefix partials**: they would double Phase 8's surface for tests-only consumers, and the cross-check fixture instantiates the right multi-index from the FD operator's tag.
- **Real-to-complex (r2c) transforms.** PocketFFT's `r2c` / `c2r` for forward and backward 3D transforms on `Field<double>`. Output spectrum has shape `(Nx/2 + 1) × Ny × Nz` (Hermitian-symmetric, x-axis half-spectrum on the i-fastest axis — matches Phase 1's locked-in storage layout). ~2× cheaper memory and time vs c2c. The y- and z-axes remain full `Ny` / `Nz` modes (no Hermitian compression on non-rfft axes).
- **Wavenumber convention.** $k_x[n_x] = 2\pi n_x / (N_x h_x)$ for $n_x \in [0, N_x/2]$ (rfft half-spectrum, positive frequencies only — the negative half is the Hermitian conjugate). For y- and z-axes (full spectra): the standard `numpy.fft`-style sign convention, $k_a[n_a] = 2\pi n_a^{\text{signed}} / (N_a h_a)$ where $n_a^{\text{signed}} = n_a$ if $n_a < N_a/2$ else $n_a - N_a$. Lives in a new header `spectral/Wavenumbers.hpp` that the operator headers consume.
- **Nyquist-mode zeroing.** For each axis $a$ with **odd** per-axis derivative count $N_a$, the Nyquist mode (index $n_a = N_a/2$) is set to zero before applying $(ik_a)^{N_a}$. Even-order axes keep the Nyquist mode untouched. Rationale: the centered-difference Nyquist mode of an anti-symmetric (odd-rank) derivative is identically zero in the continuum but quantized to $\pm \pi / h$ on the grid, producing a sign ambiguity that aliases as a high-frequency artifact under inverse-transform. Zeroing it is the standard convention for odd-rank spectral derivatives. The choice is encoded once in the `spectral::partial` helper, not duplicated per operator.
- **CI cross-check fixture covers every Phase 1–8 FD operator.** Coverage list:
  - Phase 1: `diff::laplacian<2>`, `<4>` (2 cases).
  - Phase 2: `diff::gradient<{2,4}>`, `diff::divergence<{2,4}>` (4 cases).
  - Phase 8: every `dxx`, `dyy`, `dzz`, `dxy`, `dxz`, `dyz`, all 10 rank-3 partials, all 15 rank-4 partials, `biharmonic` × `Order ∈ {2, 4}` = 64 cases.
  - **Total: 70 cross-check tests.** Each runs the convergence sweep `N ∈ {16, 32, 64, 128}` on a smooth trigonometric manufactured solution (a single mode well below the Nyquist on every grid in the sweep) and asserts log-log slope of FD–FFT discrepancy in $[\text{Order} - 0.5, \text{Order} + 0.5]$. The fixture is bundled by Phase-and-rank macro families so the test code stays under ~250 LOC despite the large coverage.
  - Becomes a **permanent CI gate** from Phase 9 forward (per roadmap).
- **Sweep grid sizes: `N ∈ {16, 32, 64, 128}`.** User pick. FFT at $N = 128$ is ~32 MB / transform; the most expensive cross-check is the order-4 rank-4 partial at $N = 128$. Estimated total CI time for the 70 cross-check tests: ~90–120 s.
- **Build flag `-DGRIDCALC_USE_FFT=ON` (default on).** Off → `spectral/` headers are excluded from compilation, the cross-check fixture is `#ifdef`-guarded out, and the cross-check CI gate is skipped. The FD code path is unaffected. The flag is propagated to the `gridcalc` INTERFACE target via `target_compile_definitions(gridcalc INTERFACE GRIDCALC_USE_FFT)` when on.
- **Explicit smooth-mode manufactured solution.** $\psi(\mathbf{x}) = \sin(\mathbf{k}\cdot\mathbf{x})$ with $\mathbf{k} = (1, 1, 1)$ on the periodic $[0, 2\pi)^3$ box for the cross-check sweeps. The mode is well below the Nyquist on every grid in the sweep ($k_a = 1 \ll N_a/2$ even at $N = 16$), so the spectral path has zero aliasing error and the FD–FFT discrepancy is dominated by the FD truncation error.
- **Version bump: `0.9.0`.** New public namespace `gridcalc::spectral`; new build option; vendored dependency. SemVer minor bump.

## Non-goals

- **The spectral path is not a primary differentiation engine.** Per `mission.md` and the roadmap: FD operators stay the production code. `gridcalc::spectral` exists for verification only. Documented in headers and in the User Guide note.
- **No spectral counterparts on non-orthogonal Bravais grids.** Phase 15.
- **No spectral counterparts on multi-atom-basis grids.** Phase 16.
- **No 31 named `spectral::*` d-prefix partials.** The generic `spectral::partial<Nx, Ny, Nz>(field)` covers them.
- **No first-class FFT operations exposed beyond `Fft.hpp`.** Inverse transforms, padding, zero-extension, real-input convenience APIs not used by the cross-check are out of scope.
- **No spectral solver path.** Phase 19's implicit-scheme spec round will weigh "FFT-diagonalization solver" again; the lightweight verification-only PocketFFT vendor here does not commit the project to that direction.

## Deferred questions

- **PocketFFT bump cadence.** Pin first; revisit if a CVE or correctness fix lands upstream.
- **Spectral counterparts on non-orthogonal / multi-atom-basis grids.** Phases 15 and 16.
- **Persistent FFT plan caching.** Phase 20 (performance pass) — PocketFFT's plan creation is fast enough at Phase 9's grid sizes that creating a plan per call is fine.
- **`fp32` / `complex<float>` spectral path.** Phase 20 if benchmarks justify; the production FD path is `double` only at this point.
