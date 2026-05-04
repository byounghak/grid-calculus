// Cahn–Hilliard coarsening demo.
//
// Solves the conservative phase-field equation
//
//     d_t psi = M * lap( -psi + psi^3 - kappa * lap(psi) )
//
// on a periodic 3D box with explicit RK4 (Phase 6 solve::integrate) and
// the standard 2nd-order finite-difference Laplacian (Phase 1
// diff::laplacian). The chemical potential mu = -psi + psi^3 - kappa lap psi
// is the variational derivative of the standard Landau-double-well +
// (kappa/2)|grad psi|^2 free energy. See the User Guide chapter 12 and
// Developer Note chapter 11 for derivations and references.
//
// Run with --help to list CLI options. VTK STRUCTURED_POINTS snapshots
// are written every --snapshot-every steps to --out-dir; per-snapshot
// diagnostics (free energy F, gradient energy E_grad, mean <psi>,
// max|psi|) are printed to stdout.

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/solve/Integrator.hpp>

#include <common/cahn_hilliard.hpp>
#include <common/vtk_writer.hpp>

namespace {

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Vec3d;
namespace ch = gridcalc::examples::cahn_hilliard;

struct Options {
    int Nx = 64;
    double dt = 0.01;
    int n_steps = 10000;
    int snapshot_every = 500;
    double kappa = 1.0;
    double mobility = 1.0;
    std::uint64_t seed = 42;
    std::string out_dir = "./out";
};

void printHelp() {
    std::cout << "Usage: cahn_hilliard [options]\n"
              << "  --n-x N              grid extent on each axis (cubic). Default 64.\n"
              << "  --dt DT              time step. Default 0.01.\n"
              << "  --n-steps N          total RK4 steps. Default 10000.\n"
              << "  --snapshot-every N   write VTK every N steps. Default 500.\n"
              << "  --kappa K            gradient-energy coefficient. Default 1.0.\n"
              << "  --mobility M         CH mobility. Default 1.0.\n"
              << "  --seed S             RNG seed for the IC. Default 42.\n"
              << "  --out-dir DIR        output directory. Default ./out.\n"
              << "  --help               print this message and exit.\n";
}

double parseDouble(const std::string& flag, const std::string& value) {
    try {
        return std::stod(value);
    } catch (const std::exception&) {
        throw std::invalid_argument("cahn_hilliard: " + flag + " expects a numeric argument; got '" +
                                    value + "'");
    }
}

int parseInt(const std::string& flag, const std::string& value) {
    try {
        return std::stoi(value);
    } catch (const std::exception&) {
        throw std::invalid_argument("cahn_hilliard: " + flag + " expects an integer; got '" +
                                    value + "'");
    }
}

Options parseOptions(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        std::string flag = argv[i];
        if (flag == "--help") {
            printHelp();
            std::exit(0);
        }
        if (i + 1 >= argc) {
            throw std::invalid_argument("cahn_hilliard: " + flag + " requires a value");
        }
        std::string val = argv[++i];
        if (flag == "--n-x") o.Nx = parseInt(flag, val);
        else if (flag == "--dt") o.dt = parseDouble(flag, val);
        else if (flag == "--n-steps") o.n_steps = parseInt(flag, val);
        else if (flag == "--snapshot-every") o.snapshot_every = parseInt(flag, val);
        else if (flag == "--kappa") o.kappa = parseDouble(flag, val);
        else if (flag == "--mobility") o.mobility = parseDouble(flag, val);
        else if (flag == "--seed") o.seed = static_cast<std::uint64_t>(std::stoull(val));
        else if (flag == "--out-dir") o.out_dir = val;
        else throw std::invalid_argument("cahn_hilliard: unknown option '" + flag + "'");
    }
    if (o.Nx <= 0 || o.dt <= 0.0 || o.n_steps < 0 || o.snapshot_every <= 0 ||
        o.kappa <= 0.0 || o.mobility <= 0.0) {
        throw std::invalid_argument("cahn_hilliard: --n-x, --dt, --kappa, --mobility must be > 0; "
                                    "--n-steps must be >= 0; --snapshot-every must be > 0");
    }
    return o;
}

// Uniform random IC on (-amplitude, amplitude), then mean-subtracted so
// the discrete sum over psi is zero to round-off. CH dynamics conserves
// the spatial mean exactly on a periodic domain, so the simulation
// stays at <psi> = 0.
Field<double> makeRandomIc(const Grid& grid, std::uint64_t seed, double amplitude) {
    Field<double> psi(grid);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> dist(-amplitude, amplitude);
    const std::size_t n = grid.getSize();
    double sum = 0.0;
    for (std::size_t idx = 0; idx < n; ++idx) {
        const double v = dist(rng);
        psi.data()[idx] = v;
        sum += v;
    }
    const double mean = sum / static_cast<double>(n);
    for (std::size_t idx = 0; idx < n; ++idx) {
        psi.data()[idx] -= mean;
    }
    return psi;
}

double computeMaxAbs(const Field<double>& psi) {
    double m = 0.0;
    const std::size_t n = psi.getSize();
    const double* p = psi.data();
    for (std::size_t idx = 0; idx < n; ++idx) {
        const double v = std::abs(p[idx]);
        if (v > m) m = v;
    }
    return m;
}

void writeSnapshot(const Field<double>& psi, const std::string& out_dir, int step) {
    std::ostringstream filename;
    filename << "psi_" << std::setw(6) << std::setfill('0') << step << ".vtk";
    const auto path = std::filesystem::path(out_dir) / filename.str();
    gridcalc::examples::writeVtkStructuredPoints(psi, path.string(), "psi");
}

void printDiagnostics(int step, double t, const Field<double>& psi, double kappa) {
    const double F = ch::computeFreeEnergy(psi, kappa);
    const double E_grad = ch::computeGradientEnergy(psi);
    const double mean = ch::computeMean(psi);
    const double max_abs = computeMaxAbs(psi);
    std::printf("step=%6d  t=%9.4f  F=%+.6e  E_grad=%.6e  <psi>=%+.3e  max|psi|=%.4f\n",
                step, t, F, E_grad, mean, max_abs);
    std::fflush(stdout);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options opts = parseOptions(argc, argv);

        std::filesystem::create_directories(opts.out_dir);

        Grid grid(opts.Nx, opts.Nx, opts.Nx, Vec3d(1.0, 1.0, 1.0));
        Field<double> psi = makeRandomIc(grid, opts.seed, 0.05);

        std::cout << "Cahn-Hilliard demo  N=" << opts.Nx << "^3"
                  << "  dt=" << opts.dt
                  << "  n_steps=" << opts.n_steps
                  << "  snapshot_every=" << opts.snapshot_every
                  << "  kappa=" << opts.kappa
                  << "  M=" << opts.mobility
                  << "  seed=" << opts.seed
                  << "  out_dir=" << opts.out_dir
                  << '\n';

        printDiagnostics(0, 0.0, psi, opts.kappa);
        writeSnapshot(psi, opts.out_dir, 0);

        const auto rhs = [&](const Field<double>& y) {
            return ch::computeRhs(y, opts.mobility, opts.kappa);
        };

        const auto t_start = std::chrono::steady_clock::now();
        int step = 0;
        while (step < opts.n_steps) {
            const int chunk = std::min(opts.snapshot_every, opts.n_steps - step);
            gridcalc::solve::integrate(psi, rhs, opts.dt, chunk, gridcalc::solve::RK4{});
            step += chunk;
            const double t = static_cast<double>(step) * opts.dt;
            printDiagnostics(step, t, psi, opts.kappa);
            writeSnapshot(psi, opts.out_dir, step);
        }
        const auto t_end = std::chrono::steady_clock::now();
        const double wall = std::chrono::duration<double>(t_end - t_start).count();
        std::cout << "wall_clock=" << wall << "s\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "cahn_hilliard: " << e.what() << '\n';
        return 1;
    }
}
