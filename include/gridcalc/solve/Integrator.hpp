/// \file
/// \brief Aggregator header for the time-integrator family.
///
/// Includes both `ExplicitEuler.hpp` and `RK4.hpp`. Each per-integrator
/// header defines its own tag struct (e.g., `gridcalc::solve::ExplicitEuler`,
/// `gridcalc::solve::RK4`) and the corresponding `solve::integrate(...)`
/// overload, so callers can use either dispatcher after including just
/// this one header.
///
/// `solve::integrate(psi, rhs, dt, n_steps, Tag{})` is the canonical
/// generic entry point; the existing free function `solve::explicitEuler`
/// (Phase 5) is preserved as a thin wrapper for backward compatibility.
/// \since 0.6.0

#pragma once

#include <gridcalc/solve/ExplicitEuler.hpp>
#include <gridcalc/solve/RK4.hpp>
