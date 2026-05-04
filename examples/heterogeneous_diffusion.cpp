// Heterogeneous-D diffusion demo (Phase 14).
//
// Solves the conservative-form heat equation
//
//     d_t psi = div(D(x) grad psi)
//
// on a periodic 3D box with Phase 14's FVM cell-flux discretization
// (`fvm::cellLaplacian`) and Phase 6's tag-dispatched
// `solve::integrate(...)`. The diffusivity `D(x) = D0 + eps cos(2 pi
// x / L)` is smooth and strictly positive; the IC is an isotropic
// Gaussian peak. Watch the Gaussian spread anisotropically through
// the heterogeneous medium: faster in high-D regions, slower in
// low-D regions.
//
// Run with --help for the CLI option list. VTK STRUCTURED_POINTS
// snapshots are written every --snapshot-every steps to --out-dir
// (reuses Phase 12's `examples/common/vtk_writer.hpp`).

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/solve/Diffusion.hpp>

#include <common/heterogeneous_diffusion.hpp>
#include <common/vtk_writer.hpp>

namespace {

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
namespace hd = gridcalc::examples::heterogeneous_diffusion;

constexpr double kPi = 3.14159265358979323846;

struct Options {
    int Nx = 64;
    double dt = 0.0;  // 0 = "auto from CFL"
    int n_steps = 2000;
    int snapshot_every = 100;
    double D0 = 0.10;
    double D_amplitude = 0.05;
    double gaussian_sigma = 0.6;
    char D_axis = 'x';
    std::string integrator = "rk4";
    std::string out_dir = "./out";
};

void printHelp() {
    std::cout
        << "Usage: heterogeneous_diffusion [options]\n"
        << "  --n-x N              cubic grid extent. Default 64.\n"
        << "  --dt DT              explicit time step (0 = auto from CFL). Default 0.\n"
        << "  --n-steps N          total integrator steps. Default 2000.\n"
        << "  --snapshot-every N   VTK snapshot every N steps. Default 100.\n"
        << "  --d-mean D0          mean diffusivity. Default 0.10.\n"
        << "  --d-amplitude EPS    cosine-amplitude of D, must be < d-mean. Default 0.05.\n"
        << "  --d-axis x|y|z       axis along which D varies. Default x.\n"
        << "  --gaussian-sigma S   Gaussian-IC standard deviation. Default 0.6.\n"
        << "  --integrator euler|rk4   time integrator. Default rk4.\n"
        << "  --out-dir DIR        VTK output directory. Default ./out.\n"
        << "  --help               print this message and exit.\n";
}

double parseDouble(const std::string& flag, const std::string& v) {
    try { return std::stod(v); } catch (...) {
        throw std::invalid_argument("heterogeneous_diffusion: " + flag +
                                    " expects a numeric argument; got '" + v + "'");
    }
}

int parseInt(const std::string& flag, const std::string& v) {
    try { return std::stoi(v); } catch (...) {
        throw std::invalid_argument("heterogeneous_diffusion: " + flag +
                                    " expects an integer; got '" + v + "'");
    }
}

Options parseOptions(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        std::string flag = argv[i];
        if (flag == "--help") { printHelp(); std::exit(0); }
        if (i + 1 >= argc) {
            throw std::invalid_argument("heterogeneous_diffusion: " + flag + " requires a value");
        }
        std::string val = argv[++i];
        if (flag == "--n-x") o.Nx = parseInt(flag, val);
        else if (flag == "--dt") o.dt = parseDouble(flag, val);
        else if (flag == "--n-steps") o.n_steps = parseInt(flag, val);
        else if (flag == "--snapshot-every") o.snapshot_every = parseInt(flag, val);
        else if (flag == "--d-mean") o.D0 = parseDouble(flag, val);
        else if (flag == "--d-amplitude") o.D_amplitude = parseDouble(flag, val);
        else if (flag == "--d-axis") {
            if (val != "x" && val != "y" && val != "z") {
                throw std::invalid_argument("heterogeneous_diffusion: --d-axis must be x|y|z");
            }
            o.D_axis = val[0];
        }
        else if (flag == "--gaussian-sigma") o.gaussian_sigma = parseDouble(flag, val);
        else if (flag == "--integrator") {
            if (val != "euler" && val != "rk4") {
                throw std::invalid_argument("heterogeneous_diffusion: --integrator must be euler|rk4");
            }
            o.integrator = val;
        }
        else if (flag == "--out-dir") o.out_dir = val;
        else throw std::invalid_argument("heterogeneous_diffusion: unknown option '" + flag + "'");
    }
    if (o.Nx <= 0 || o.n_steps < 0 || o.snapshot_every <= 0 ||
        o.D0 <= 0.0 || o.gaussian_sigma <= 0.0 ||
        std::abs(o.D_amplitude) >= o.D0) {
        throw std::invalid_argument("heterogeneous_diffusion: invalid CLI options "
                                    "(check Nx > 0, n_steps >= 0, D0 > 0, |amp| < D0)");
    }
    return o;
}

void writeSnapshot(const Field<double>& psi, const std::string& out_dir, int step) {
    std::ostringstream filename;
    filename << "psi_" << std::setw(6) << std::setfill('0') << step << ".vtk";
    const auto path = std::filesystem::path(out_dir) / filename.str();
    gridcalc::examples::writeVtkStructuredPoints(psi, path.string(), "psi");
}

void printDiagnostics(int step, double t, const Field<double>& psi) {
    const double mass = hd::computeMass(psi);
    const double l2sq = hd::computeL2Squared(psi);
    const double psi_min = hd::computeMin(psi);
    const double psi_max = hd::computeMax(psi);
    std::printf("step=%6d  t=%9.4f  mass=%+.6e  L2^2=%+.6e  min=%+.3e  max=%+.4f\n",
                step, t, mass, l2sq, psi_min, psi_max);
    std::fflush(stdout);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options opts = parseOptions(argc, argv);
        std::filesystem::create_directories(opts.out_dir);

        const double L = 2.0 * kPi;
        const double h = L / opts.Nx;
        Grid grid(opts.Nx, opts.Nx, opts.Nx, Vec3d(h, h, h));

        const Vec3d center(L / 2.0, L / 2.0, L / 2.0);
        Field<double> psi = hd::makeGaussianIc(grid, opts.gaussian_sigma, center);
        Field<double> D = hd::makeHeterogeneousD(grid, opts.D0, opts.D_amplitude, opts.D_axis);

        // Auto-pick dt from the CFL bound if the user passed --dt 0.
        const double D_max = opts.D0 + std::abs(opts.D_amplitude);
        const double cfl_bound =
            (opts.integrator == "rk4" ? 0.6963 : 0.5) * h * h / (3.0 * D_max);
        const double dt = (opts.dt > 0.0) ? opts.dt : 0.4 * cfl_bound;

        std::cout << "Heterogeneous-D diffusion demo  N=" << opts.Nx << "^3"
                  << "  dt=" << dt
                  << "  cfl_bound=" << cfl_bound
                  << "  n_steps=" << opts.n_steps
                  << "  D0=" << opts.D0
                  << "  D_amp=" << opts.D_amplitude
                  << "  D_axis=" << opts.D_axis
                  << "  sigma=" << opts.gaussian_sigma
                  << "  integrator=" << opts.integrator
                  << "  out=" << opts.out_dir << '\n';

        printDiagnostics(0, 0.0, psi);
        writeSnapshot(psi, opts.out_dir, 0);

        const auto t_start = std::chrono::steady_clock::now();
        int step = 0;
        while (step < opts.n_steps) {
            const int chunk = std::min(opts.snapshot_every, opts.n_steps - step);
            if (opts.integrator == "rk4") {
                gridcalc::solve::diffuse(psi, D, dt, chunk, gridcalc::solve::RK4{});
            } else {
                gridcalc::solve::diffuse(psi, D, dt, chunk, gridcalc::solve::ExplicitEuler{});
            }
            step += chunk;
            const double t = static_cast<double>(step) * dt;
            printDiagnostics(step, t, psi);
            writeSnapshot(psi, opts.out_dir, step);
        }
        const double wall =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
        std::cout << "wall_clock=" << wall << "s\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "heterogeneous_diffusion: " << e.what() << '\n';
        return 1;
    }
}
