// Linear elastic-energy demo (Phase 15).
//
// Computes the linear-elastic energy
//
//     F = integral( 0.5 * sigma : epsilon, dV )
//
// of a periodic uniaxial dilation displacement
//
//     u_axis(x) = alpha * sin(k0 * x_axis),   k0 = 1
//
// on the standard `[0, 2 pi]^3` periodic box, using Phase 15's
// `func::evaluate(Field<Vec3d>, ...)` and the Phase 13 rank-2 gradient
// `M(i, j) = d_j u_i`. The closed-form reference
//
//     F_ref = 2 pi^3 (lambda + 2 mu) alpha^2
//
// is printed alongside the discrete value so the user can see the
// O(h^2) error scale with `--n-x`.
//
// Run with --help for the CLI option list. A single VTK snapshot of
// the axial strain `epsilon_axis,axis` is written to `--out-dir` for
// inspection in ParaView / VisIt.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/diff/Gradient.hpp>
#include <gridcalc/func/Functional.hpp>

#include <common/elastic_energy.hpp>
#include <common/vtk_writer.hpp>

namespace {

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Mat3d;
using gridcalc::core::Vec3d;
namespace ee = gridcalc::examples::elastic_energy;
namespace func = gridcalc::func;

constexpr double kPi = 3.14159265358979323846;

struct Options {
    int Nx = 64;
    // Material parameters --- two routes:
    //   (a) supply E and nu (default --young 1, --nu 0.3 -> Lamé 0.5769, 0.3846)
    //   (b) override either Lamé constant directly via --lambda / --mu
    bool lambda_set = false;
    bool mu_set = false;
    double E = 1.0;
    double nu = 0.3;
    double lambda = 0.0;
    double mu = 0.0;
    double alpha = 0.01;
    int axis = 0;
    std::string out_dir = "./out";
};

void printHelp() {
    std::cout
        << "Usage: elastic_energy [options]\n"
        << "  --n-x N           cubic grid extent. Default 64.\n"
        << "  --young E         Young's modulus (used with --nu). Default 1.0.\n"
        << "  --nu NU           Poisson's ratio (used with --young). Default 0.3.\n"
        << "  --lambda L        First Lamé constant; overrides --young/--nu route.\n"
        << "  --mu M            Second Lamé constant (shear modulus); overrides --young/--nu.\n"
        << "  --alpha A         Peak strain amplitude of the dilation wave. Default 0.01.\n"
        << "  --axis x|y|z      Direction of the dilation wave. Default x.\n"
        << "  --out-dir DIR     Output directory for the VTK strain snapshot. Default ./out.\n"
        << "  --help            Print this message and exit.\n";
}

double parseDouble(const std::string& flag, const std::string& v) {
    try { return std::stod(v); } catch (...) {
        throw std::invalid_argument("elastic_energy: " + flag +
                                    " expects a numeric argument; got '" + v + "'");
    }
}

int parseInt(const std::string& flag, const std::string& v) {
    try { return std::stoi(v); } catch (...) {
        throw std::invalid_argument("elastic_energy: " + flag +
                                    " expects an integer; got '" + v + "'");
    }
}

int parseAxis(const std::string& v) {
    if (v == "x") return 0;
    if (v == "y") return 1;
    if (v == "z") return 2;
    throw std::invalid_argument("elastic_energy: --axis must be x|y|z (got '" + v + "')");
}

Options parseOptions(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        std::string flag = argv[i];
        if (flag == "--help") { printHelp(); std::exit(0); }
        if (i + 1 >= argc) {
            throw std::invalid_argument("elastic_energy: " + flag + " requires a value");
        }
        std::string val = argv[++i];
        if (flag == "--n-x") o.Nx = parseInt(flag, val);
        else if (flag == "--young") o.E = parseDouble(flag, val);
        else if (flag == "--nu") o.nu = parseDouble(flag, val);
        else if (flag == "--lambda") { o.lambda = parseDouble(flag, val); o.lambda_set = true; }
        else if (flag == "--mu") { o.mu = parseDouble(flag, val); o.mu_set = true; }
        else if (flag == "--alpha") o.alpha = parseDouble(flag, val);
        else if (flag == "--axis") o.axis = parseAxis(val);
        else if (flag == "--out-dir") o.out_dir = val;
        else throw std::invalid_argument("elastic_energy: unknown option '" + flag + "'");
    }
    if (o.Nx <= 0) {
        throw std::invalid_argument("elastic_energy: --n-x must be positive");
    }
    if (!(o.lambda_set && o.mu_set)) {
        const auto lame = ee::lameFromYoungPoisson(o.E, o.nu);
        if (!o.lambda_set) o.lambda = lame.first;
        if (!o.mu_set) o.mu = lame.second;
    }
    if (!(o.lambda + 2.0 * o.mu > 0.0) || !(o.mu > 0.0)) {
        throw std::invalid_argument(
            "elastic_energy: Lamé constants must satisfy mu > 0 and lambda + 2 mu > 0");
    }
    return o;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options opts = parseOptions(argc, argv);
        std::filesystem::create_directories(opts.out_dir);

        const double L = 2.0 * kPi;
        const double h = L / opts.Nx;
        Grid grid(opts.Nx, opts.Nx, opts.Nx, Vec3d(h, h, h));

        const Field<Vec3d> u =
            ee::makeUniaxialPeriodicDisplacement(grid, opts.alpha, opts.axis);

        const auto t_start = std::chrono::steady_clock::now();
        const double F_h = func::evaluate(u, [&opts](const Vec3d&, const Mat3d& grad_u) {
            return ee::linearElasticEnergyDensity(opts.lambda, opts.mu, grad_u);
        });
        const double F_ref =
            ee::analyticalLinearElasticEnergy(grid, opts.alpha, opts.lambda, opts.mu, opts.axis);
        const auto t_end = std::chrono::steady_clock::now();

        const double rel = std::abs(F_h - F_ref) / std::abs(F_ref);
        const double seconds =
            std::chrono::duration<double>(t_end - t_start).count();

        // Snapshot of the axial strain epsilon_{axis,axis} = d u_axis / d x_axis.
        const Field<Mat3d> grad_u = gridcalc::diff::gradient(u);
        Field<double> eps_axial(grid);
        for (int k = 0; k < grid.getNz(); ++k) {
            for (int j = 0; j < grid.getNy(); ++j) {
                for (int i = 0; i < grid.getNx(); ++i) {
                    eps_axial(i, j, k) = grad_u(i, j, k)(opts.axis, opts.axis);
                }
            }
        }
        const auto path = std::filesystem::path(opts.out_dir) / "elastic_energy.vtk";
        gridcalc::examples::writeVtkStructuredPoints(
            eps_axial, path.string(), "epsilon_axial");

        std::printf("Elastic-energy demo (Phase 15):\n");
        std::printf("  N = %d   h = %.6f   axis = %c\n", opts.Nx, h,
                    "xyz"[opts.axis]);
        std::printf("  lambda = %.6f   mu = %.6f   alpha = %.6f\n",
                    opts.lambda, opts.mu, opts.alpha);
        std::printf("  F_h     = %+.10e\n", F_h);
        std::printf("  F_ref   = %+.10e   (closed form, 2 pi^3 (lambda + 2 mu) alpha^2)\n",
                    F_ref);
        std::printf("  rel_err = %.3e   wall = %.3f s\n", rel, seconds);
        std::printf("  snapshot: %s\n", path.string().c_str());
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "elastic_energy: error -- %s\n", e.what());
        return 1;
    }
}
