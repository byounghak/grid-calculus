// Phase 15 follow-up #1, audit fix #1: prevent dangling Leaf /
// IdentityField from rvalue Field<T> / Grid wraps.
//
// Before 0.15.1, `expr::field(diff::gradient(u))` wrapped the
// temporary Field<Mat3d> returned by `gradient`, destroyed it at the
// end of the full-expression statement, and left the Leaf with a
// dangling raw pointer. Subsequent `materialize` / `func::integrate`
// walks dereferenced freed memory.
//
// All checks below are *compile-time* `static_assert`s --- the
// failure mode is now a `= delete`d overload, surfacing at compile
// time. The single runtime `TEST` body exists so `gtest_discover_tests`
// registers a name for the file (and so a regression that accidentally
// re-enables the rvalue path would still produce a failing build, not
// a silently-skipped test file).

#include <type_traits>
#include <utility>

#include <gtest/gtest.h>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>
#include <gridcalc/tensor/Expressions.hpp>

namespace {

using gridcalc::core::Field;
using gridcalc::core::Grid;
using gridcalc::core::Mat3d;
using gridcalc::core::Vec3d;
namespace expr = gridcalc::tensor::expr;

// SFINAE detector for `expr::field(declval<Arg>())`. Deleted overloads
// participate in overload resolution; selecting one causes the call
// expression to be ill-formed, which `void_t`-style SFINAE detects in
// the immediate context. Defined as a struct (not a lambda) so it
// stays C++17.
template <class T, class Arg, class = void>
struct FieldFactoryCallable : std::false_type {};

template <class T, class Arg>
struct FieldFactoryCallable<
    T, Arg, std::void_t<decltype(expr::field(std::declval<Arg>()))>>
    : std::true_type {};

template <class T, class Arg>
inline constexpr bool kFieldFactoryCallable = FieldFactoryCallable<T, Arg>::value;

template <class Arg, class = void>
struct IdentityFieldFactoryCallable : std::false_type {};

template <class Arg>
struct IdentityFieldFactoryCallable<
    Arg, std::void_t<decltype(expr::identityField(std::declval<Arg>()))>>
    : std::true_type {};

template <class Arg>
inline constexpr bool kIdentityFieldFactoryCallable =
    IdentityFieldFactoryCallable<Arg>::value;

// --- Constructor matrix --------------------------------------------------
// Lvalue references are accepted; rvalues and prvalues are rejected.

static_assert(std::is_constructible_v<expr::Leaf<double>, const Field<double>&>,
              "Leaf<double> must accept const Field<double>&");
static_assert(!std::is_constructible_v<expr::Leaf<double>, Field<double>&&>,
              "Leaf<double> must reject Field<double>&&");
static_assert(!std::is_constructible_v<expr::Leaf<double>, Field<double>>,
              "Leaf<double> must reject prvalue Field<double>");

static_assert(std::is_constructible_v<expr::Leaf<Vec3d>, const Field<Vec3d>&>,
              "Leaf<Vec3d> must accept const Field<Vec3d>&");
static_assert(!std::is_constructible_v<expr::Leaf<Vec3d>, Field<Vec3d>&&>,
              "Leaf<Vec3d> must reject Field<Vec3d>&&");
static_assert(!std::is_constructible_v<expr::Leaf<Vec3d>, Field<Vec3d>>,
              "Leaf<Vec3d> must reject prvalue Field<Vec3d>");

static_assert(std::is_constructible_v<expr::Leaf<Mat3d>, const Field<Mat3d>&>,
              "Leaf<Mat3d> must accept const Field<Mat3d>&");
static_assert(!std::is_constructible_v<expr::Leaf<Mat3d>, Field<Mat3d>&&>,
              "Leaf<Mat3d> must reject Field<Mat3d>&&");
static_assert(!std::is_constructible_v<expr::Leaf<Mat3d>, Field<Mat3d>>,
              "Leaf<Mat3d> must reject prvalue Field<Mat3d>");

static_assert(std::is_constructible_v<expr::IdentityField, const Grid&>,
              "IdentityField must accept const Grid&");
static_assert(!std::is_constructible_v<expr::IdentityField, Grid&&>,
              "IdentityField must reject Grid&&");
static_assert(!std::is_constructible_v<expr::IdentityField, Grid>,
              "IdentityField must reject prvalue Grid");

// --- Factory matrix ------------------------------------------------------
// SFINAE-detected callability of the overload sets `expr::field` and
// `expr::identityField`. The lvalue overload is callable; the
// `= delete`d rvalue overload makes the call expression ill-formed.

static_assert(kFieldFactoryCallable<Mat3d, const Field<Mat3d>&>,
              "expr::field<Mat3d> must accept const Field<Mat3d>&");
static_assert(!kFieldFactoryCallable<Mat3d, Field<Mat3d>&&>,
              "expr::field<Mat3d> must reject Field<Mat3d>&&");

static_assert(kFieldFactoryCallable<Vec3d, const Field<Vec3d>&>,
              "expr::field<Vec3d> must accept const Field<Vec3d>&");
static_assert(!kFieldFactoryCallable<Vec3d, Field<Vec3d>&&>,
              "expr::field<Vec3d> must reject Field<Vec3d>&&");

static_assert(kFieldFactoryCallable<double, const Field<double>&>,
              "expr::field<double> must accept const Field<double>&");
static_assert(!kFieldFactoryCallable<double, Field<double>&&>,
              "expr::field<double> must reject Field<double>&&");

static_assert(kIdentityFieldFactoryCallable<const Grid&>,
              "expr::identityField must accept const Grid&");
static_assert(!kIdentityFieldFactoryCallable<Grid&&>,
              "expr::identityField must reject Grid&&");

// --- Runtime sanity ------------------------------------------------------

constexpr double kPi = 3.14159265358979323846;

TEST(TensorLeafLifetimeTest, LvalueFactoriesStillWork) {
    const double h = 2.0 * kPi / 4;
    const Grid g(4, 4, 4, Vec3d(h, h, h));
    const Field<Mat3d> M(g, Mat3d::Identity());
    const auto leaf = expr::field(M);
    const auto idf = expr::identityField(g);
    // `Field` holds its grid by value; `leaf.grid()` references the
    // field's internal copy. Compare extents and cell sizes instead.
    EXPECT_EQ(leaf.grid().getNx(), g.getNx());
    EXPECT_EQ(leaf.grid().getNy(), g.getNy());
    EXPECT_EQ(leaf.grid().getNz(), g.getNz());
    EXPECT_EQ(&idf.grid(), &g);  // IdentityField stores the original Grid by reference.
    EXPECT_TRUE(leaf.evalAt(0, 0, 0).isApprox(Mat3d::Identity(), 1e-15));
    EXPECT_TRUE(idf.evalAt(0, 0, 0).isApprox(Mat3d::Identity(), 1e-15));
}

}  // namespace
