/// \file
/// \brief Project-local Eigen typedefs shared across `gridcalc` headers.
/// \since 0.2.0

#pragma once

#include <Eigen/Core>

namespace gridcalc::core {

/// \brief Project-local alias for the 3-component Cartesian vector type.
///
/// Originally declared in `core/Grid.hpp` (Phase 1) and relocated here in
/// Phase 2 once `diff/Gradient.hpp` and `diff/Divergence.hpp` also needed it.
/// The fully-qualified name `gridcalc::core::Vec3d` is unchanged from
/// Phase 1.
/// \since 0.1.0
using Vec3d = Eigen::Vector3d;

}  // namespace gridcalc::core
