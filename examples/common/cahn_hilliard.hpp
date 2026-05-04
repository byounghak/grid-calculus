// Cahn–Hilliard RHS and energy diagnostics for the demo and its
// acceptance test. The free energy form is the standard
// Landau-double-well plus gradient term:
//
//     F[psi] = \int [ -psi^2/2 + psi^4/4 + (kappa/2) |grad psi|^2 ] dV ,
//
// whose variational derivative (the chemical potential) is
//
//     mu = delta F / delta psi = -psi + psi^3 - kappa lap psi .
//
// The conservative dynamics
//
//     d_t psi = M lap (mu)
//
// is the canonical model B / Cahn–Hilliard equation for binary
// phase-separation. See `docs/developer-note/chapters/11-cahn-hilliard.tex`
// for the derivation, stability analysis, and references.

#pragma once

#include <cstddef>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/diff/Laplacian.hpp>
#include <gridcalc/func/Functional.hpp>

namespace gridcalc::examples::cahn_hilliard {

/// Evaluates the Cahn–Hilliard right-hand side
/// `M * lap(-psi + psi^3 - kappa * lap(psi))`. Two `diff::laplacian`
/// calls per evaluation, one allocation each. Returned `Field` is the
/// time derivative `d_t psi` ready to be plugged into
/// `solve::integrate(... RK4{})`.
inline core::Field<double> computeRhs(const core::Field<double>& psi,
                                      double M,
                                      double kappa) {
    auto lap_psi = diff::laplacian(psi);

    core::Field<double> mu(psi.getGrid());
    const std::size_t n = psi.getSize();
    const double* p_psi = psi.data();
    const double* p_lap = lap_psi.data();
    double* p_mu = mu.data();
    for (std::size_t idx = 0; idx < n; ++idx) {
        const double p = p_psi[idx];
        p_mu[idx] = -p + p * p * p - kappa * p_lap[idx];
    }

    auto rhs = diff::laplacian(mu);
    double* p_rhs = rhs.data();
    for (std::size_t idx = 0; idx < n; ++idx) {
        p_rhs[idx] *= M;
    }
    return rhs;
}

/// Total free energy `F[psi] = \int [ -psi^2/2 + psi^4/4 + (kappa/2)|grad psi|^2 ] dV`.
/// Implemented via `func::evaluate` to exercise the Phase 4 / Phase 11
/// 2-arg arity directly.
inline double computeFreeEnergy(const core::Field<double>& psi, double kappa) {
    return func::evaluate(psi, [kappa](double p, const core::Vec3d& g) {
        return -0.5 * p * p + 0.25 * p * p * p * p + 0.5 * kappa * g.squaredNorm();
    });
}

/// Gradient (interfacial) energy `\int |grad psi|^2 dV`. Decreases
/// monotonically once the system enters the coarsening regime; serves
/// as the primary statistic for the acceptance test.
inline double computeGradientEnergy(const core::Field<double>& psi) {
    return func::evaluate(psi, [](double, const core::Vec3d& g) {
        return g.squaredNorm();
    });
}

/// Mean field `<psi> = (1/V) \int psi dV`. CH dynamics conserves it
/// exactly in the continuum and to round-off in the discretized
/// periodic implementation; the acceptance test asserts this.
inline double computeMean(const core::Field<double>& psi) {
    const auto& grid = psi.getGrid();
    const double total_volume = static_cast<double>(grid.getSize()) * grid.getCellVolume();
    const double integral = func::evaluate(psi, [](double p) { return p; });
    return integral / total_volume;
}

}  // namespace gridcalc::examples::cahn_hilliard
