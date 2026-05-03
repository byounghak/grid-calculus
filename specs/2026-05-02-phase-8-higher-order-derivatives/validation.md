# Validation: Phase 8 — Higher-order spatial derivatives (3rd, 4th)

## Build

- [ ] CI green on Ubuntu GCC and Ubuntu Clang.
- [ ] Warnings clean under `-Wall -Wextra -Wpedantic -Wconversion`.
- [ ] `clang-format` and `clang-tidy` pass (no new regressions vs the codebase-wide hygiene gap noted in earlier STATUS entries).
- [ ] `gridcalc` target remains INTERFACE — no new source files.
- [ ] `ThirdDerivative<3>`, `FourthDerivative<6>`, and any other unsupported `Order` continue to produce a compile error if instantiated (primary templates stay undefined; verified by inspection — no positive test exists because a successful compile of an unsupported instantiation would itself be the failure).

## Tests

- [ ] **Weight-table audits** — `EXPECT_DOUBLE_EQ` passes for every weight in `ThirdDerivative<2>`, `ThirdDerivative<4>`, `FourthDerivative<2>`, `FourthDerivative<4>`.
- [ ] **Convergence sweeps** — for every `(operator, Order)` pair across the 31 named partials and `biharmonic`, the log-log slope on `N ∈ {16, 32, 64}` lies in `[Order − 0.5, Order + 0.5]`. Total: **64 sweep tests**.
  - Rank-2 (12 tests): `dxx`, `dyy`, `dzz`, `dxy`, `dxz`, `dyz` × `Order ∈ {2, 4}`.
  - Rank-3 (20 tests): `d3dx3`, `d3dy3`, `d3dz3`, `d3dx2dy`, `d3dxdy2`, `d3dx2dz`, `d3dxdz2`, `d3dy2dz`, `d3dydz2`, `d3dxdydz` × `{2, 4}`.
  - Rank-4 (30 tests): `d4dx4`, `d4dy4`, `d4dz4`, `d4dx3dy`, `d4dxdy3`, `d4dx3dz`, `d4dxdz3`, `d4dy3dz`, `d4dydz3`, `d4dx2dy2`, `d4dx2dz2`, `d4dy2dz2`, `d4dx2dydz`, `d4dxdy2dz`, `d4dxdydz2` × `{2, 4}`.
  - Biharmonic (2 tests): `biharmonic` × `{2, 4}`.
- [ ] **Roadmap acceptance — biharmonic recovers $|\mathbf{k}|^4$.** `BiharmonicTest.RecoversKToTheFourthOrder2` and `BiharmonicTest.RecoversKToTheFourthOrder4`: for $\psi = \sin(\mathbf{k}\cdot\mathbf{x})$ on $[0, 2\pi)^3$ with at least two non-trivial $\mathbf{k}$ vectors, `biharmonic<Order>(psi)` recovers $|\mathbf{k}|^4 \psi$ pointwise within the order-`Order` truncation error at `N = 64`.
- [ ] **`d4` alias parity.** `d4<Order>(psi)` returns a `Field<double>` that is exactly equal (`EXPECT_EQ` on every grid point) to `biharmonic<Order>(psi)` for the same input. Two cases (`Order ∈ {2, 4}`).
- [ ] **All Phase 1–7 tests still pass** without modification.

## Documentation

- [ ] Every new public symbol carries a Doxygen block including `\brief`, `\param`, `\returns`, `\tparam Order`, `\since 0.8.0`. Covers: 31 named partial-derivative functions, `biharmonic`, `d4` alias, `stencil::ThirdDerivative<2>`, `stencil::ThirdDerivative<4>`, `stencil::FourthDerivative<2>`, `stencil::FourthDerivative<4>`, the `detail::mixedDerivative` helper.
- [ ] User Guide note added at `docs/user-guide/notes/phase-8-higher-order-derivatives.md` covering: the full naming convention, the 31 named functions plus `biharmonic`/`d4`, calling syntax with `Order`, a worked $\nabla^4$ example with the spectral expectation cross-check.
- [ ] Developer Note added at `docs/developer-note/notes/phase-8-higher-order-derivatives.md` with the five-section structure (Theory · Math derivation · Algorithm · Design decisions · References) per `specs/CLAUDE.md` step 4d. References section non-empty: at least one external citation (DOI / textbook bibliographic entry / permanent URL).
- [ ] `CHANGELOG.md` entry under `0.8.0` lists the new stencil tables (`ThirdDerivative<Order>`, `FourthDerivative<Order>`), the 31 named partial-derivative functions, `biharmonic` + `d4` alias, and the convergence-test coverage.
- [ ] `STATUS.md` updated: Phase 8 row marked Done; Phase 9 stated as the new Next action. Add d-prefix-naming-convention bullet and Fornberg-deferral bullet to "Decisions worth knowing."
- [ ] CI `render-docs` job produces fresh aggregate PDFs that include the Phase 8 chapter in both books.

## Performance

- N/A at Phase 8. The most stencil-intensive operator at order 4 is `d4dx2dy2<4>` with a 5×5 = 25-point cross stencil (per `(i, j, k)` the helper reads up to $\prod_a (2 r_a + 1)$ points; for an order-4 rank-(2, 2, 0) op that is 5 × 5 × 1 = 25 reads). Expected ~25× the order-2 single-axis operator cost; documented in the Developer Note. Phase 20 owns optimization.
