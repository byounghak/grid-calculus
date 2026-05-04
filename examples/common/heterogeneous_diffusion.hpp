// IC + diagnostics helpers shared between the heterogeneous-D
// diffusion demo (`examples/heterogeneous_diffusion.cpp`) and its CI
// acceptance gate (`test/heterogeneous_diffusion_test.cpp`).
//
// Mirrors the Phase 12 `examples/common/cahn_hilliard.hpp` pattern: the
// helpers live under `examples/` (not `include/gridcalc/`) because they
// are not part of the public library API; the test target's
// `target_include_directories` entry from Phase 12 makes them
// reachable from `test/` as `#include <common/...>`.

#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>

namespace gridcalc::examples::heterogeneous_diffusion {

/// Builds an isotropic Gaussian peak `psi(x) = exp(-|x - center|^2 /
/// (2 sigma^2))` on the given grid. The Gaussian is *not* normalized;
/// peak amplitude is exactly 1.0 at `center`. Boundary effects are
/// negligible when `sigma << L_box / 4` and the Gaussian is centered
/// at least ~3 sigma from any edge.
inline core::Field<double> makeGaussianIc(const core::Grid& grid,
                                          double sigma,
                                          const core::Vec3d& center) {
    return core::Field<double>(grid, [&](double x, double y, double z) {
        const double dx = x - center(0);
        const double dy = y - center(1);
        const double dz = z - center(2);
        return std::exp(-(dx * dx + dy * dy + dz * dz) / (2.0 * sigma * sigma));
    });
}

/// Builds a smooth heterogeneous-D field `D(x) = D0 + eps * cos(2 pi x
/// / L)` along the chosen axis (`'x' | 'y' | 'z'`). Caller must pass
/// `eps < D0` so the result is strictly positive everywhere — the
/// `solve::diffuse(psi, D, ...)` heterogeneous-D overload requires
/// `D > 0` at every cell. Throws `std::invalid_argument` on the
/// positivity violation.
inline core::Field<double> makeHeterogeneousD(const core::Grid& grid,
                                              double D0,
                                              double eps,
                                              char axis) {
    if (D0 <= 0.0) {
        throw std::invalid_argument(
            "makeHeterogeneousD: D0 must be > 0; got " + std::to_string(D0));
    }
    if (std::abs(eps) >= D0) {
        throw std::invalid_argument(
            "makeHeterogeneousD: |eps| must be < D0 to keep D > 0; got eps=" +
            std::to_string(eps) + ", D0=" + std::to_string(D0));
    }

    const double L = (axis == 'x') ? grid.getNx() * grid.getCellSize()(0)
                   : (axis == 'y') ? grid.getNy() * grid.getCellSize()(1)
                                   : grid.getNz() * grid.getCellSize()(2);
    const double k = 2.0 * 3.14159265358979323846 / L;

    return core::Field<double>(grid, [=](double x, double y, double z) {
        const double s = (axis == 'x') ? x : (axis == 'y') ? y : z;
        return D0 + eps * std::cos(k * s);
    });
}

/// Total mass `Σψ * cell_volume`. CH dynamics conserves it on a
/// periodic domain; the demo prints this every snapshot to confirm.
inline double computeMass(const core::Field<double>& psi) {
    double sum = 0.0;
    for (double v : psi) sum += v;
    return sum * psi.getGrid().getCellVolume();
}

/// L²-norm-squared `∫ψ² dV`. Strictly decreasing under heterogeneous
/// diffusion with `D > 0` (the textbook entropy-like Lyapunov
/// statement: `d/dt ∫ψ²/2 dV = -∫D|∇ψ|² dV ≤ 0`). The demo prints this
/// every snapshot; the CI gate asserts it decreases between chunks.
inline double computeL2Squared(const core::Field<double>& psi) {
    double sum = 0.0;
    for (double v : psi) sum += v * v;
    return sum * psi.getGrid().getCellVolume();
}

/// Pointwise minimum of `psi`. Used by the demo + acceptance test
/// to verify positivity preservation: a positive Gaussian IC stays
/// `>= 0` to round-off under heterogeneous-D diffusion with `D > 0`.
inline double computeMin(const core::Field<double>& psi) {
    double m = psi.data()[0];
    for (double v : psi) {
        if (v < m) m = v;
    }
    return m;
}

/// Pointwise maximum of `psi`. Should monotonically decrease as the
/// Gaussian spreads (max-principle / dissipativity).
inline double computeMax(const core::Field<double>& psi) {
    double m = psi.data()[0];
    for (double v : psi) {
        if (v > m) m = v;
    }
    return m;
}

}  // namespace gridcalc::examples::heterogeneous_diffusion
