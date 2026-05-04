// Hand-written legacy ASCII VTK STRUCTURED_POINTS writer for `Field<double>`.
//
// Writes a single snapshot to a `.vtk` file readable by ParaView and VisIt
// without any new build dependency. The format is documented at
// https://kitware.github.io/vtk-examples/site/VTKFileFormats/ ; we emit the
// "Simple Legacy Formats" variant with `ASCII` data encoding.
//
// This helper lives under examples/ because it is not part of the public
// gridcalc library API. Production users wanting binary or XML VTK should
// integrate vtkIO directly.

#pragma once

#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>

#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>

namespace gridcalc::examples {

/// Write `field` as a legacy ASCII VTK STRUCTURED_POINTS file.
///
/// File layout: the standard VTK 3.0 header, followed by `DIMENSIONS Nx Ny
/// Nz`, `ORIGIN 0 0 0`, `SPACING hx hy hz`, then a `POINT_DATA` block with
/// one `SCALARS <field_name> double 1` array. Values are emitted in the
/// `Field`'s i-fastest storage order (x varies fastest) — exactly the order
/// VTK expects for STRUCTURED_POINTS, so no transposition is needed.
///
/// On open failure the function throws `std::runtime_error`. The caller is
/// responsible for ensuring the parent directory exists.
inline void writeVtkStructuredPoints(const core::Field<double>& field,
                                     const std::string& path,
                                     const std::string& field_name) {
    const core::Grid& grid = field.getGrid();
    const auto& h = grid.getCellSize();

    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("writeVtkStructuredPoints: cannot open '" + path +
                                 "' for writing");
    }
    out.precision(17);

    out << "# vtk DataFile Version 3.0\n";
    out << "gridcalc field: " << field_name << "\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << grid.getNx() << ' ' << grid.getNy() << ' ' << grid.getNz() << '\n';
    out << "ORIGIN 0 0 0\n";
    out << "SPACING " << h(0) << ' ' << h(1) << ' ' << h(2) << '\n';
    out << "POINT_DATA " << grid.getSize() << '\n';
    out << "SCALARS " << field_name << " double 1\n";
    out << "LOOKUP_TABLE default\n";

    const std::size_t n = grid.getSize();
    const double* data = field.data();
    for (std::size_t idx = 0; idx < n; ++idx) {
        out << data[idx] << '\n';
    }

    if (!out) {
        throw std::runtime_error("writeVtkStructuredPoints: write to '" + path + "' failed");
    }
}

}  // namespace gridcalc::examples
