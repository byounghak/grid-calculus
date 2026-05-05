/// \file
/// \brief Contraction expression-template (ET) layer for `Field<T>` arithmetic
///        with `T \in \{double, core::Vec3d, core::Mat3d\}`.
///
/// Phase 13 shipped the eager primitives `tensor::trace`,
/// `tensor::singleContract`, and `tensor::doubleContract`, each of which
/// allocates a fresh `Field<...>` for its output. Phase 15 adds a parallel
/// lazy surface under `tensor::expr::*`: every operator builds an
/// expression node that exposes `evalAt(i, j, k)` and propagates a `Grid`
/// reference up the tree. `tensor::expr::materialize(E)` walks all cells
/// once and stores the result into a fresh `Field<value_type<E>>`;
/// `func::integrate(E, tag)` consumes a `double`-valued expression and
/// reduces in a single pass without materializing rank-2 intermediates.
///
/// The motivating workload is the linear elastic energy
/// $F = \int \tfrac12 \boldsymbol{\sigma}:\boldsymbol{\varepsilon}\,dV$:
/// the ET layer evaluates the per-point Lamé closed form
/// $\tfrac12 (\lambda \mathrm{tr}(\boldsymbol{\varepsilon})^2 + 2\mu\,\boldsymbol{\varepsilon}:\boldsymbol{\varepsilon})$
/// in a fused loop with no `Field<core::Mat3d>` intermediates ever
/// allocated for $\boldsymbol{\sigma}$, $\boldsymbol{\varepsilon}$, or any
/// linear combination thereof.
///
/// Phase 13's eager `tensor::*` primitives are not modified --- the new
/// nodes are a parallel surface, not a replacement.
/// \since 0.15.0

#pragma once

#include <cstddef>
#include <type_traits>
#include <utility>

#include <gridcalc/core/EigenAliases.hpp>
#include <gridcalc/core/Field.hpp>
#include <gridcalc/core/Grid.hpp>

namespace gridcalc::tensor::expr {

/// \brief Trait registering a type as a contraction-ET node.
///
/// Specialized to `std::true_type` for every node type defined below
/// (`Leaf`, `IdentityField`, `Scaled`, `Plus`, `Minus`, `Negate`, `Trace`,
/// `Sym`, `SingleContract`, `DoubleContract`). Used by the operator
/// overloads at the bottom of this header to gate ET-node arithmetic
/// without colliding with arithmetic on unrelated types.
/// \since 0.15.0
template <class T>
struct is_expr : std::false_type {};

/// \brief Convenience variable template for `is_expr<T>::value`.
/// \since 0.15.0
template <class T>
inline constexpr bool is_expr_v = is_expr<std::remove_cv_t<std::remove_reference_t<T>>>::value;

// ---------------------------------------------------------------------------
//  Leaf nodes
// ---------------------------------------------------------------------------

/// \brief Lazy view of a `core::Field<T>` for `T \in \{double, core::Vec3d, core::Mat3d\}`.
///
/// Holds a non-owning pointer to the wrapped field; the field must outlive
/// any expression containing this leaf. Created via the `expr::field`
/// factory.
/// \tparam T  Element type of the wrapped `Field`.
/// \since 0.15.0
template <class T>
class Leaf {
 public:
  using value_type = T;

  /// \brief Wraps a reference to the given field.
  /// \param f  The field to view. Must outlive `*this`.
  /// \since 0.15.0
  explicit Leaf(const core::Field<T>& f) noexcept : _field(&f) {}

  /// \brief Returns the grid the underlying field lives on.
  /// \since 0.15.0
  const core::Grid& grid() const noexcept { return _field->getGrid(); }

  /// \brief Returns the wrapped field's value at `(i, j, k)`.
  /// \since 0.15.0
  const T& evalAt(int i, int j, int k) const noexcept { return (*_field)(i, j, k); }

 private:
  const core::Field<T>* _field;
};

template <class T>
struct is_expr<Leaf<T>> : std::true_type {};

/// \brief Wraps a `core::Field<T>` as a `Leaf<T>` expression node.
/// \param f  Field to wrap (lifetime contract: must outlive any expression containing the leaf).
/// \returns A `Leaf<T>` referencing `f`.
/// \since 0.15.0
template <class T>
inline Leaf<T> field(const core::Field<T>& f) noexcept {
  return Leaf<T>{f};
}

/// \brief Lazy broadcast of `core::Mat3d::Identity()` over every grid cell.
///
/// Carries no per-cell storage --- `evalAt` returns the identity matrix
/// directly. Used to inject the Lamé `\lambda \mathrm{tr}(\varepsilon) \mathbf{I}`
/// term into stress expressions without materializing a constant
/// `Field<core::Mat3d>`.
/// \since 0.15.0
class IdentityField {
 public:
  using value_type = core::Mat3d;

  /// \brief Constructs a broadcast identity over the given grid.
  /// \param g  Grid the resulting expression is defined on. Must outlive `*this`.
  /// \since 0.15.0
  explicit IdentityField(const core::Grid& g) noexcept : _grid(&g) {}

  /// \brief Returns the grid the broadcast lives on.
  /// \since 0.15.0
  const core::Grid& grid() const noexcept { return *_grid; }

  /// \brief Returns `core::Mat3d::Identity()` regardless of `(i, j, k)`.
  /// \since 0.15.0
  core::Mat3d evalAt(int /*i*/, int /*j*/, int /*k*/) const noexcept {
    return core::Mat3d::Identity();
  }

 private:
  const core::Grid* _grid;
};

template <>
struct is_expr<IdentityField> : std::true_type {};

/// \brief Returns a broadcast `Mat3d::Identity()` expression over the given grid.
/// \param g  Grid the resulting expression is defined on (referenced; must outlive returned node).
/// \returns An `IdentityField` expression.
/// \since 0.15.0
inline IdentityField identityField(const core::Grid& g) noexcept {
  return IdentityField{g};
}

// ---------------------------------------------------------------------------
//  Combinator nodes
// ---------------------------------------------------------------------------

/// \brief Pointwise scalar multiplication `c * E`.
/// \tparam E  Inner expression type.
/// \since 0.15.0
template <class E>
class Scaled {
 public:
  using value_type = typename E::value_type;

  /// \brief Constructs `c * inner` (stored by value).
  /// \since 0.15.0
  Scaled(double c, E inner) noexcept(std::is_nothrow_move_constructible_v<E>)
      : _c(c), _inner(std::move(inner)) {}

  /// \brief Returns the inner expression's grid.
  /// \since 0.15.0
  const core::Grid& grid() const noexcept { return _inner.grid(); }

  /// \brief Returns `c * inner.evalAt(i, j, k)`.
  /// \since 0.15.0
  value_type evalAt(int i, int j, int k) const noexcept {
    return _c * _inner.evalAt(i, j, k);
  }

 private:
  double _c;
  E _inner;
};

template <class E>
struct is_expr<Scaled<E>> : std::true_type {};

/// \brief Pointwise sum of two same-`value_type` expressions.
/// \since 0.15.0
template <class L, class R>
class Plus {
 public:
  static_assert(std::is_same_v<typename L::value_type, typename R::value_type>,
                "tensor::expr::Plus: operands must have the same value_type");
  using value_type = typename L::value_type;

  /// \brief Constructs `lhs + rhs` (operands stored by value).
  /// \since 0.15.0
  Plus(L lhs, R rhs) noexcept(std::is_nothrow_move_constructible_v<L> &&
                              std::is_nothrow_move_constructible_v<R>)
      : _lhs(std::move(lhs)), _rhs(std::move(rhs)) {}

  /// \brief Returns the (shared) grid of the operands.
  /// \since 0.15.0
  const core::Grid& grid() const noexcept { return _lhs.grid(); }

  /// \brief Returns `lhs.evalAt(...) + rhs.evalAt(...)`.
  /// \since 0.15.0
  value_type evalAt(int i, int j, int k) const noexcept {
    return _lhs.evalAt(i, j, k) + _rhs.evalAt(i, j, k);
  }

 private:
  L _lhs;
  R _rhs;
};

template <class L, class R>
struct is_expr<Plus<L, R>> : std::true_type {};

/// \brief Pointwise difference of two same-`value_type` expressions.
/// \since 0.15.0
template <class L, class R>
class Minus {
 public:
  static_assert(std::is_same_v<typename L::value_type, typename R::value_type>,
                "tensor::expr::Minus: operands must have the same value_type");
  using value_type = typename L::value_type;

  /// \brief Constructs `lhs - rhs` (operands stored by value).
  /// \since 0.15.0
  Minus(L lhs, R rhs) noexcept(std::is_nothrow_move_constructible_v<L> &&
                               std::is_nothrow_move_constructible_v<R>)
      : _lhs(std::move(lhs)), _rhs(std::move(rhs)) {}

  /// \brief Returns the (shared) grid of the operands.
  /// \since 0.15.0
  const core::Grid& grid() const noexcept { return _lhs.grid(); }

  /// \brief Returns `lhs.evalAt(...) - rhs.evalAt(...)`.
  /// \since 0.15.0
  value_type evalAt(int i, int j, int k) const noexcept {
    return _lhs.evalAt(i, j, k) - _rhs.evalAt(i, j, k);
  }

 private:
  L _lhs;
  R _rhs;
};

template <class L, class R>
struct is_expr<Minus<L, R>> : std::true_type {};

// ---------------------------------------------------------------------------
//  Tensor primitives (rank-2 only --- rank-3+ deferred per Phase 15 scope)
// ---------------------------------------------------------------------------

/// \brief Pointwise trace of a `core::Mat3d`-valued expression.
///
/// `Trace(M)(i, j, k) = M.evalAt(i, j, k).trace() = M(0, 0) + M(1, 1) + M(2, 2)`.
/// Under the Phase 13 rank-2 gradient convention `M(i, j) = \partial_j v_i`,
/// the trace gives the divergence of the underlying vector field.
/// \since 0.15.0
template <class E>
class Trace {
 public:
  static_assert(std::is_same_v<typename E::value_type, core::Mat3d>,
                "tensor::expr::Trace requires a Mat3d-valued expression");
  using value_type = double;

  /// \brief Constructs `trace(inner)` (inner stored by value).
  /// \since 0.15.0
  explicit Trace(E inner) noexcept(std::is_nothrow_move_constructible_v<E>)
      : _inner(std::move(inner)) {}

  /// \brief Returns the inner expression's grid.
  /// \since 0.15.0
  const core::Grid& grid() const noexcept { return _inner.grid(); }

  /// \brief Returns `inner.evalAt(i, j, k).trace()`.
  /// \since 0.15.0
  double evalAt(int i, int j, int k) const noexcept {
    return _inner.evalAt(i, j, k).trace();
  }

 private:
  E _inner;
};

template <class E>
struct is_expr<Trace<E>> : std::true_type {};

/// \brief Returns a lazy `Trace` of a `Mat3d`-valued expression.
/// \param e  Expression with `value_type == core::Mat3d`.
/// \returns A `Trace<std::decay_t<E>>` expression node.
/// \since 0.15.0
template <class E, std::enable_if_t<is_expr_v<E>, int> = 0>
inline auto trace(E&& e) {
  return Trace<std::decay_t<E>>{std::forward<E>(e)};
}

/// \brief Pointwise symmetric part `\tfrac12 (M + M^\top)` of a `Mat3d`-valued expression.
/// \since 0.15.0
template <class E>
class Sym {
 public:
  static_assert(std::is_same_v<typename E::value_type, core::Mat3d>,
                "tensor::expr::Sym requires a Mat3d-valued expression");
  using value_type = core::Mat3d;

  /// \brief Constructs `sym(inner)`.
  /// \since 0.15.0
  explicit Sym(E inner) noexcept(std::is_nothrow_move_constructible_v<E>)
      : _inner(std::move(inner)) {}

  /// \brief Returns the inner expression's grid.
  /// \since 0.15.0
  const core::Grid& grid() const noexcept { return _inner.grid(); }

  /// \brief Returns `0.5 * (inner.evalAt(...) + inner.evalAt(...).transpose())`.
  /// \since 0.15.0
  core::Mat3d evalAt(int i, int j, int k) const noexcept {
    const core::Mat3d m = _inner.evalAt(i, j, k);
    return 0.5 * (m + m.transpose());
  }

 private:
  E _inner;
};

template <class E>
struct is_expr<Sym<E>> : std::true_type {};

/// \brief Returns a lazy `Sym` of a `Mat3d`-valued expression.
/// \param e  Expression with `value_type == core::Mat3d`.
/// \returns A `Sym<std::decay_t<E>>` expression node.
/// \since 0.15.0
template <class E, std::enable_if_t<is_expr_v<E>, int> = 0>
inline auto sym(E&& e) {
  return Sym<std::decay_t<E>>{std::forward<E>(e)};
}

/// \brief Pointwise single (matrix-product) contraction of two `Mat3d`-valued expressions.
///
/// `(A \cdot B)(i, j) = A(i, k) B(k, j)` (sum over `k`); identical to the
/// Phase 13 eager `tensor::singleContract`.
/// \since 0.15.0
template <class L, class R>
class SingleContract {
 public:
  static_assert(std::is_same_v<typename L::value_type, core::Mat3d>,
                "tensor::expr::SingleContract: lhs must have value_type == Mat3d");
  static_assert(std::is_same_v<typename R::value_type, core::Mat3d>,
                "tensor::expr::SingleContract: rhs must have value_type == Mat3d");
  using value_type = core::Mat3d;

  /// \brief Constructs `singleContract(lhs, rhs)`.
  /// \since 0.15.0
  SingleContract(L lhs, R rhs) noexcept(std::is_nothrow_move_constructible_v<L> &&
                                        std::is_nothrow_move_constructible_v<R>)
      : _lhs(std::move(lhs)), _rhs(std::move(rhs)) {}

  /// \brief Returns the (shared) grid of the operands.
  /// \since 0.15.0
  const core::Grid& grid() const noexcept { return _lhs.grid(); }

  /// \brief Returns `lhs.evalAt(...) * rhs.evalAt(...)` (Eigen `Matrix3d::operator*`).
  /// \since 0.15.0
  core::Mat3d evalAt(int i, int j, int k) const noexcept {
    return _lhs.evalAt(i, j, k) * _rhs.evalAt(i, j, k);
  }

 private:
  L _lhs;
  R _rhs;
};

template <class L, class R>
struct is_expr<SingleContract<L, R>> : std::true_type {};

/// \brief Returns a lazy single-contraction (matrix product) of two `Mat3d`-valued expressions.
/// \param a  Left operand.
/// \param b  Right operand.
/// \since 0.15.0
template <class L, class R, std::enable_if_t<is_expr_v<L> && is_expr_v<R>, int> = 0>
inline auto singleContract(L&& a, R&& b) {
  return SingleContract<std::decay_t<L>, std::decay_t<R>>{std::forward<L>(a),
                                                          std::forward<R>(b)};
}

/// \brief Pointwise full (double) contraction of two `Mat3d`-valued expressions.
///
/// `A : B = A(i, j) B(i, j)` (sum over both indices). Identical to the
/// Phase 13 eager `tensor::doubleContract`.
/// \since 0.15.0
template <class L, class R>
class DoubleContract {
 public:
  static_assert(std::is_same_v<typename L::value_type, core::Mat3d>,
                "tensor::expr::DoubleContract: lhs must have value_type == Mat3d");
  static_assert(std::is_same_v<typename R::value_type, core::Mat3d>,
                "tensor::expr::DoubleContract: rhs must have value_type == Mat3d");
  using value_type = double;

  /// \brief Constructs `doubleContract(lhs, rhs)`.
  /// \since 0.15.0
  DoubleContract(L lhs, R rhs) noexcept(std::is_nothrow_move_constructible_v<L> &&
                                        std::is_nothrow_move_constructible_v<R>)
      : _lhs(std::move(lhs)), _rhs(std::move(rhs)) {}

  /// \brief Returns the (shared) grid of the operands.
  /// \since 0.15.0
  const core::Grid& grid() const noexcept { return _lhs.grid(); }

  /// \brief Returns `(lhs.evalAt(...).array() * rhs.evalAt(...).array()).sum()`.
  /// \since 0.15.0
  double evalAt(int i, int j, int k) const noexcept {
    const core::Mat3d a = _lhs.evalAt(i, j, k);
    const core::Mat3d b = _rhs.evalAt(i, j, k);
    return (a.array() * b.array()).sum();
  }

 private:
  L _lhs;
  R _rhs;
};

template <class L, class R>
struct is_expr<DoubleContract<L, R>> : std::true_type {};

/// \brief Returns a lazy full-contraction `A : B` of two `Mat3d`-valued expressions.
/// \param a  Left operand.
/// \param b  Right operand.
/// \since 0.15.0
template <class L, class R, std::enable_if_t<is_expr_v<L> && is_expr_v<R>, int> = 0>
inline auto doubleContract(L&& a, R&& b) {
  return DoubleContract<std::decay_t<L>, std::decay_t<R>>{std::forward<L>(a),
                                                          std::forward<R>(b)};
}

// ---------------------------------------------------------------------------
//  Operator overloads (gated on `is_expr_v`)
// ---------------------------------------------------------------------------

/// \brief Unary minus of an expression: returns `Scaled<E>(-1.0, e)`.
/// \since 0.15.0
template <class E, std::enable_if_t<is_expr_v<E>, int> = 0>
inline auto operator-(E&& e) {
  return Scaled<std::decay_t<E>>{-1.0, std::forward<E>(e)};
}

/// \brief Left scalar multiplication `c * E`.
/// \since 0.15.0
template <class E, std::enable_if_t<is_expr_v<E>, int> = 0>
inline auto operator*(double c, E&& e) {
  return Scaled<std::decay_t<E>>{c, std::forward<E>(e)};
}

/// \brief Right scalar multiplication `E * c`.
/// \since 0.15.0
template <class E, std::enable_if_t<is_expr_v<E>, int> = 0>
inline auto operator*(E&& e, double c) {
  return Scaled<std::decay_t<E>>{c, std::forward<E>(e)};
}

/// \brief Right scalar division `E / c`.
/// \since 0.15.0
template <class E, std::enable_if_t<is_expr_v<E>, int> = 0>
inline auto operator/(E&& e, double c) {
  return Scaled<std::decay_t<E>>{1.0 / c, std::forward<E>(e)};
}

/// \brief Pointwise sum of two expressions (must share `value_type`).
/// \since 0.15.0
template <class L, class R, std::enable_if_t<is_expr_v<L> && is_expr_v<R>, int> = 0>
inline auto operator+(L&& l, R&& r) {
  return Plus<std::decay_t<L>, std::decay_t<R>>{std::forward<L>(l), std::forward<R>(r)};
}

/// \brief Pointwise difference of two expressions (must share `value_type`).
/// \since 0.15.0
template <class L, class R, std::enable_if_t<is_expr_v<L> && is_expr_v<R>, int> = 0>
inline auto operator-(L&& l, R&& r) {
  return Minus<std::decay_t<L>, std::decay_t<R>>{std::forward<L>(l), std::forward<R>(r)};
}

// ---------------------------------------------------------------------------
//  Materialization
// ---------------------------------------------------------------------------

/// \brief Walks every grid cell, evaluates `expr.evalAt(...)`, and stores the result
///        into a fresh `core::Field<value_type<E>>`.
///
/// The escape hatch when a caller needs the materialized field as
/// downstream input. For reductions (`func::integrate(expr, tag)`),
/// prefer the fused-loop overload in `gridcalc/func/Integrate.hpp`,
/// which avoids the temporary allocation when the result is being
/// summed away anyway.
/// \tparam E    Expression-node type (`is_expr_v<E>` must be `true`).
/// \param expr  The expression to materialize.
/// \returns A new `core::Field<value_type<E>>` carrying the per-cell values.
/// \since 0.15.0
template <class E, std::enable_if_t<is_expr_v<E>, int> = 0>
inline core::Field<typename std::decay_t<E>::value_type> materialize(const E& expr) {
  using T = typename std::decay_t<E>::value_type;
  const core::Grid& g = expr.grid();
  core::Field<T> out(g);
  const int Nx = g.getNx();
  const int Ny = g.getNy();
  const int Nz = g.getNz();
  for (int k = 0; k < Nz; ++k) {
    for (int j = 0; j < Ny; ++j) {
      for (int i = 0; i < Nx; ++i) {
        out(i, j, k) = expr.evalAt(i, j, k);
      }
    }
  }
  return out;
}

}  // namespace gridcalc::tensor::expr
