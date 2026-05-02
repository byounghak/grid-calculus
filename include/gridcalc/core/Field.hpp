/// \file
/// \brief Declares `gridcalc::core::Field`, a 3D scalar (or arbitrary `T`)
///        field stored on a `Grid`.
/// \since 0.1.0

#pragma once

#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

#include <gridcalc/core/Grid.hpp>
#include <gridcalc/core/IndexPolicy.hpp>

namespace gridcalc::core {

/// \brief A 3D field of values of type `T` stored on a `Grid` with a chosen
///        out-of-range index policy.
///
/// Storage layout is i-fastest:
/// `linear_index = i + Nx * (j + Ny * k)`. The `Grid` is held by value so a
/// `Field` is self-contained (no lifetime contract on an external `Grid`).
///
/// `Policy` defaults to `IndexPolicy::Periodic`. `operator()(i, j, k)` runs
/// `Policy::wrap` on each axis before indexing, so callers may pass negative
/// or `>= N` coordinates and receive the policy-defined value (the wrapped
/// neighbor for `Periodic`). Phase 18 introduces `Dirichlet` and `Neumann`
/// policies.
///
/// \tparam T       Element type. `double` is the only instantiation used in
///                 Phase 1; non-scalar `T` (e.g. `Vec3d`) lands in Phase 13.
/// \tparam Policy  Stateless out-of-range index policy. Must expose
///                 `static constexpr int wrap(int i, int N) noexcept`.
/// \since 0.1.0
template <class T, class Policy = IndexPolicy::Periodic>
class Field {
 public:
  /// \brief Element type alias for STL-style introspection.
  /// \since 0.1.0
  using value_type = T;
  /// \brief Policy alias, e.g. for selecting compatible operators.
  /// \since 0.1.0
  using policy_type = Policy;

  /// \brief Constructs a zero-initialized `Field` on the given `Grid`.
  /// \param grid  The sampling grid. Stored by value.
  /// \since 0.1.0
  explicit Field(const Grid& grid) : _grid(grid), _data(grid.getSize()) {}

  /// \brief Constructs a `Field` filled with a single broadcast `value`.
  /// \param grid   The sampling grid. Stored by value.
  /// \param value  Value broadcast to every grid point.
  /// \since 0.1.0
  Field(const Grid& grid, const T& value) : _grid(grid), _data(grid.getSize(), value) {}

  /// \brief Constructs a `Field` by sampling a Cartesian function `f(x, y, z)`.
  ///
  /// Sample point `(i, j, k)` is evaluated at `(i*hx, j*hy, k*hz)`. The grid
  /// origin is `(0, 0, 0)`; no half-cell offset is applied. Disabled via
  /// SFINAE when `f` is not invocable as `f(double, double, double) -> T` so
  /// that `Field(grid, value)` keeps winning for value-typed second args.
  /// \param grid  The sampling grid.
  /// \param f     Callable `(double, double, double) -> T` (or convertible).
  /// \since 0.1.0
  template <class F,
            std::enable_if_t<std::is_invocable_r_v<T, F&, double, double, double>, int> = 0>
  Field(const Grid& grid, F&& f) : _grid(grid), _data(grid.getSize()) {
    const auto& h = _grid.getCellSize();
    for (int k = 0; k < _grid.getNz(); ++k) {
      for (int j = 0; j < _grid.getNy(); ++j) {
        for (int i = 0; i < _grid.getNx(); ++i) {
          _data[_grid.getLinearIndex(i, j, k)] = static_cast<T>(
              std::forward<F>(f)(static_cast<double>(i) * h(0),
                                 static_cast<double>(j) * h(1),
                                 static_cast<double>(k) * h(2)));
        }
      }
    }
  }

  /// \brief Mutable element access. Applies `Policy::wrap` on each axis.
  /// \param i  x-coordinate, possibly out of `[0, Nx)`.
  /// \param j  y-coordinate, possibly out of `[0, Ny)`.
  /// \param k  z-coordinate, possibly out of `[0, Nz)`.
  /// \returns Reference to the stored element at the wrapped position.
  /// \since 0.1.0
  T& operator()(int i, int j, int k) noexcept {
    return _data[_grid.getLinearIndex(Policy::wrap(i, _grid.getNx()),
                                      Policy::wrap(j, _grid.getNy()),
                                      Policy::wrap(k, _grid.getNz()))];
  }

  /// \brief Const element access. Applies `Policy::wrap` on each axis.
  /// \param i  x-coordinate, possibly out of `[0, Nx)`.
  /// \param j  y-coordinate, possibly out of `[0, Ny)`.
  /// \param k  z-coordinate, possibly out of `[0, Nz)`.
  /// \returns Const reference to the stored element at the wrapped position.
  /// \since 0.1.0
  const T& operator()(int i, int j, int k) const noexcept {
    return _data[_grid.getLinearIndex(Policy::wrap(i, _grid.getNx()),
                                      Policy::wrap(j, _grid.getNy()),
                                      Policy::wrap(k, _grid.getNz()))];
  }

  /// \brief Returns the `Grid` this `Field` was built on.
  /// \returns Const reference to the stored `Grid`.
  /// \since 0.1.0
  const Grid& getGrid() const noexcept { return _grid; }

  /// \brief Returns the total number of stored elements (`Nx * Ny * Nz`).
  /// \since 0.1.0
  std::size_t getSize() const noexcept { return _data.size(); }

  /// \brief STL-protocol shim for `getSize()`.
  /// \since 0.1.0
  std::size_t size() const noexcept { return getSize(); }

  /// \brief Returns a pointer to the underlying contiguous storage.
  /// \since 0.1.0
  T* data() noexcept { return _data.data(); }

  /// \brief Const overload of `data()`.
  /// \since 0.1.0
  const T* data() const noexcept { return _data.data(); }

  /// \brief Iterator to the first stored element.
  /// \since 0.1.0
  T* begin() noexcept { return _data.data(); }
  /// \brief Const iterator to the first stored element.
  /// \since 0.1.0
  const T* begin() const noexcept { return _data.data(); }
  /// \brief Iterator past the last stored element.
  /// \since 0.1.0
  T* end() noexcept { return _data.data() + _data.size(); }
  /// \brief Const iterator past the last stored element.
  /// \since 0.1.0
  const T* end() const noexcept { return _data.data() + _data.size(); }

 private:
  /// \brief Sampling grid. Owned by value.
  Grid _grid;
  /// \brief i-fastest contiguous storage of length `Nx * Ny * Nz`.
  std::vector<T> _data;
};

}  // namespace gridcalc::core
