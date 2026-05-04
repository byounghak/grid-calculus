#include <common/vtk_writer.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>

namespace {

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;

std::filesystem::path makeUniqueTempFile(const std::string& stem) {
    auto base = std::filesystem::temp_directory_path() / "gridcalc_vtk_writer_test";
    std::filesystem::create_directories(base);
    return base / (stem + "_" +
                   std::to_string(reinterpret_cast<std::uintptr_t>(&stem)) + ".vtk");
}

// Build a deterministic non-trivial pattern so the round-trip catches
// transposition bugs and i-fastest-vs-k-fastest mistakes.
Field<double> makePattern(int Nx, int Ny, int Nz, const Vec3d& h) {
    Grid g(Nx, Ny, Nz, h);
    Field<double> f(g);
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                f(i, j, k) = static_cast<double>(i) + 0.5 * static_cast<double>(j) +
                             0.25 * static_cast<double>(k);
            }
        }
    }
    return f;
}

TEST(VtkWriterTest, WriteThenParse_HeaderMatchesGrid) {
    const int Nx = 4, Ny = 3, Nz = 2;
    const Vec3d h(0.1, 0.2, 0.3);
    auto field = makePattern(Nx, Ny, Nz, h);

    auto path = makeUniqueTempFile("header");
    gridcalc::examples::writeVtkStructuredPoints(field, path.string(), "psi");

    std::ifstream in(path);
    ASSERT_TRUE(in.good());

    std::string line;
    std::vector<std::string> lines;
    while (std::getline(in, line)) lines.push_back(line);

    ASSERT_EQ(lines[0], "# vtk DataFile Version 3.0");
    EXPECT_EQ(lines[1], "gridcalc field: psi");
    EXPECT_EQ(lines[2], "ASCII");
    EXPECT_EQ(lines[3], "DATASET STRUCTURED_POINTS");
    EXPECT_EQ(lines[4], "DIMENSIONS 4 3 2");
    EXPECT_EQ(lines[5], "ORIGIN 0 0 0");
    EXPECT_EQ(lines[6].substr(0, 8), "SPACING ");
    EXPECT_EQ(lines[7], "POINT_DATA 24");
    EXPECT_EQ(lines[8], "SCALARS psi double 1");
    EXPECT_EQ(lines[9], "LOOKUP_TABLE default");
    EXPECT_EQ(lines.size(), static_cast<std::size_t>(10 + Nx * Ny * Nz));

    std::filesystem::remove(path);
}

TEST(VtkWriterTest, WriteThenParse_PayloadMatchesField) {
    const int Nx = 4, Ny = 3, Nz = 2;
    const Vec3d h(1.0, 1.0, 1.0);
    auto field = makePattern(Nx, Ny, Nz, h);

    auto path = makeUniqueTempFile("payload");
    gridcalc::examples::writeVtkStructuredPoints(field, path.string(), "psi");

    std::ifstream in(path);
    std::string line;
    for (int i = 0; i < 10; ++i) std::getline(in, line);  // skip header

    // Parse data values. VTK STRUCTURED_POINTS expects values in the same
    // i-fastest order our `Field` uses, so a flat read should match
    // `field.data()` element-by-element.
    std::vector<double> parsed;
    parsed.reserve(field.getSize());
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        parsed.push_back(std::stod(line));
    }

    ASSERT_EQ(parsed.size(), field.getSize());
    for (std::size_t idx = 0; idx < field.getSize(); ++idx) {
        EXPECT_DOUBLE_EQ(parsed[idx], field.data()[idx])
            << "mismatch at flat index " << idx;
    }

    std::filesystem::remove(path);
}

TEST(VtkWriterTest, WriteThenParse_FieldNameRoundTrips) {
    Grid g(2, 2, 2, Vec3d(1.0, 1.0, 1.0));
    Field<double> field(g, 0.5);

    auto path = makeUniqueTempFile("name");
    gridcalc::examples::writeVtkStructuredPoints(field, path.string(), "rho_total");

    std::ifstream in(path);
    std::string line;
    bool found_scalars_line = false;
    while (std::getline(in, line)) {
        if (line == "SCALARS rho_total double 1") {
            found_scalars_line = true;
            break;
        }
    }
    EXPECT_TRUE(found_scalars_line);

    std::filesystem::remove(path);
}

}  // namespace
