// Phase 15 elastic-energy helpers, shared between the demo binary
// (`examples/elastic_energy.cpp`) and its CI acceptance gate
// (`test/elastic_energy_test.cpp` and `test/elastic_energy_demo_test.cpp`).
//
// Mirrors the Phase 12 / Phase 14 `examples/common/*.hpp` pattern: the
// helpers live under `examples/` (not `include/gridcalc/`) because they
// are not part of the public library API; the test target's
// `target_include_directories` entry from Phase 12 makes them
// reachable from `test/` as `#include <common/...>`.

#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>

namespace gridcalc::examples::elastic_energy {

/// Returns Lamé constants `(lambda, mu)` from Young's modulus `E` and
/// Poisson's ratio `nu`. Validates `E > 0` and `nu in (-1, 0.5)` (the
/// physically admissible range for isotropic linear elasticity).
inline std::pair<double, double> lameFromYoungPoisson(double E, double nu) {
    if (!(E > 0.0)) {
        throw std::invalid_argument(
            "elastic_energy::lameFromYoungPoisson: Young's modulus E must be positive");
    }
    if (!(nu > -1.0 && nu < 0.5)) {
        throw std::invalid_argument(
            "elastic_energy::lameFromYoungPoisson: Poisson's ratio nu must lie in (-1, 1/2)");
    }
    const double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mu = E / (2.0 * (1.0 + nu));
    return {lambda, mu};
}

/// Builds the periodic uniaxial dilation displacement `u_axis(x_axis) =
/// alpha * sin(k0 * x_axis)`, with `k0 = 1` so a single sin wave fits
/// in the standard `[0, 2 pi]^3` periodic box; all other displacement
/// components are zero.
///
/// Strain (= sym(grad u)) is `diag(alpha cos(k0 x_axis), 0, 0)` permuted
/// to the chosen axis, so `trace(eps) = alpha cos(k0 x_axis)` and
/// `eps : eps = alpha^2 cos^2(k0 x_axis)`. The closed-form total elastic
/// energy on the standard box is
/// `F_ref = 2 pi^3 (lambda + 2 mu) alpha^2`, independent of `axis`.
///
/// `axis` selects 0, 1, or 2 for x, y, z. Throws `std::invalid_argument`
/// for any other value.
inline core::Field<core::Vec3d> makeUniaxialPeriodicDisplacement(const core::Grid& grid,
                                                                 double alpha,
                                                                 int axis = 0) {
    if (axis < 0 || axis > 2) {
        throw std::invalid_argument(
            "elastic_energy::makeUniaxialPeriodicDisplacement: axis must be 0, 1, or 2 (got " +
            std::to_string(axis) + ")");
    }
    constexpr double kK0 = 1.0;
    return core::Field<core::Vec3d>(grid, [alpha, axis](double x, double y, double z) {
        const double coord = (axis == 0) ? x : (axis == 1) ? y : z;
        core::Vec3d u = core::Vec3d::Zero();
        u(axis) = alpha * std::sin(kK0 * coord);
        return u;
    });
}

/// Pointwise linear-elastic energy density
/// `0.5 sigma : eps = 0.5 (lambda * trace(eps)^2 + 2 mu * eps : eps)`,
/// where `eps = sym(grad_u)`. Return value is the integrand of the
/// total elastic energy `F = integral(0.5 * sigma : eps, dV)`.
///
/// Suitable as the `func::evaluate(Field<Vec3d>, ...)` callable: pass
/// `lambda` and `mu` by capture and write the per-point body in two
/// lines. The closed-form Lamé identity sidesteps building the full
/// stress tensor `sigma = lambda * trace(eps) * I + 2 mu * eps`,
/// avoiding three Mat3d intermediates per cell.
inline double linearElasticEnergyDensity(double lambda,
                                         double mu,
                                         const core::Mat3d& grad_u) noexcept {
    const core::Mat3d eps = 0.5 * (grad_u + grad_u.transpose());
    const double tr = eps.trace();
    const double eps_double_eps = (eps.array() * eps.array()).sum();
    return 0.5 * (lambda * tr * tr + 2.0 * mu * eps_double_eps);
}

/// Closed-form total linear-elastic energy for the periodic uniaxial
/// dilation wave produced by `makeUniaxialPeriodicDisplacement`. Returns
/// `2 pi^3 (lambda + 2 mu) alpha^2` on the standard `[0, 2 pi]^3` box;
/// the value does not depend on the axis argument by symmetry. The
/// `axis` parameter is retained for explicit symmetry and so callers can
/// document which axis they used.
///
/// Throws `std::invalid_argument` if the grid is not the standard
/// `[0, 2 pi]^3` box (axis-aligned cell sizes within `1e-12` of
/// `2 pi / N_axis`); the closed form is specific to `k0 = 1` on that
/// box.
inline double analyticalLinearElasticEnergy(const core::Grid& grid,
                                            double alpha,
                                            double lambda,
                                            double mu,
                                            int axis = 0) {
    if (axis < 0 || axis > 2) {
        throw std::invalid_argument(
            "elastic_energy::analyticalLinearElasticEnergy: axis must be 0, 1, or 2 (got " +
            std::to_string(axis) + ")");
    }
    constexpr double kPi = 3.14159265358979323846;
    constexpr double kTol = 1e-12;
    const auto& h = grid.getCellSize();
    const double expected_hx = 2.0 * kPi / grid.getNx();
    const double expected_hy = 2.0 * kPi / grid.getNy();
    const double expected_hz = 2.0 * kPi / grid.getNz();
    if (std::abs(h(0) - expected_hx) > kTol || std::abs(h(1) - expected_hy) > kTol ||
        std::abs(h(2) - expected_hz) > kTol) {
        throw std::invalid_argument(
            "elastic_energy::analyticalLinearElasticEnergy: closed form requires the standard "
            "[0, 2 pi]^3 periodic box (cell sizes 2 pi / N on each axis).");
    }
    return 2.0 * kPi * kPi * kPi * (lambda + 2.0 * mu) * alpha * alpha;
}

}  // namespace gridcalc::examples::elastic_energy
