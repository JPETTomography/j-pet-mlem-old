#pragma once

#include <cmath>

/// `NYI` *Not Yet Implemented*

/// This namespace implementation is not finished, neither used across the
/// project.
namespace Geometry {
namespace {  // keeps Base private
namespace Base {

template <int N, typename F> class Vector {
 public:
  Vector(F v[N]) {
    for (int i = 0; i < N; i++) {
      v_[i] = v[i];
    }
  }

#ifndef _MSC_VER
  // MSVC does not support array initializers
  template <typename... Fs> Vector(Fs... vs) : v_{ vs... } {}
#endif

  Vector& operator+=(const Vector& rhs) {
    for (int i = 0; i < N; ++i)
      v_[i] += rhs[i];
    return *this;
  }

  Vector& operator-=(const Vector& rhs) {
    for (int i = 0; i < N; ++i)
      v_[i] -= rhs[i];
    return *this;
  }

  Vector& operator+=(F rhs) {
    for (int i = 0; i < N; ++i)
      v_[i] += rhs;
    return *this;
  }

  Vector& operator-=(F rhs) {
    for (int i = 0; i < N; ++i)
      v_[i] -= rhs;
    return *this;
  }

  Vector& operator*=(F rhs) {
    for (int i = 0; i < N; ++i)
      v_[i] *= rhs;
    return *this;
  }

  Vector& operator/=(F rhs) {
    for (int i = 0; i < N; ++i)
      v_[i] /= rhs;
    return *this;
  }

  F& operator[](int i) { return v_[i]; }
  F operator[](int i) const { return v_[i]; }

  bool operator==(const Vector& rhs) {
    for (int i = 0; i < N; ++i)
      if (rhs.v_[i] != v_[i])
        return false;
    return true;
  }

 private:
  F v_[N];
};

template <int N, typename F> class Point {
 public:
  Point(F p[N]) {
    for (int i = 0; i < N; ++i)
      p_[i] += p[i];
  }

#ifndef _MSC_VER
  // MSVC does not support array initializers
  template <typename... Fs> Point(Fs... ps) : p_{ ps... } {}
#endif

  Point& operator+=(const Vector<N, F>& rhs) {
    for (int i = 0; i < N; ++i)
      p_[i] += rhs[i];

    return *this;
  }

  Point& operator-=(const Vector<N, F>& rhs) {
    for (int i = 0; i < N; ++i)
      p_[i] -= rhs[i];
    return *this;
  }

  bool operator==(const Point& rhs) {
    for (int i = 0; i < N; ++i)
      if (rhs.p_[i] != p_[i])
        return false;
    return true;
  }

  F& operator[](int i) { return p_[i]; }
  F operator[](int i) const { return p_[i]; }

 private:
  F p_[N];
};
}  // Base
}  // private

/// \cond PRIVATE
template <int N, typename F = double> class Vector : public Base::Vector<N, F> {
 public:
  using Base::Vector<N, F>::Vector;  // inherit constructors
};
/// \endcond

template <typename F> class Vector<2, F> : public Base::Vector<2, F> {
 public:
  using Base::Vector<2, F>::Vector;  // inherit constructors

  Vector& rotate(F angle) {
    F s = std::sin(angle);
    F c = std::cos(angle);
    F x = (*this)[0];
    F y = (*this)[1];

    (*this)[0] = x * c - y * s;
    (*this)[1] = y * c + x * s;

    return *this;
  }
};

/// \cond PRIVATE
template <int N, typename F = double> class Point : public Base::Point<N, F> {
 public:
  using Base::Point<N, F>::Point;  // inherit constructors
};
/// \endcond

template <typename F> class Point<2, F> : public Base::Point<2, F> {
 public:
  using Base::Point<2, F>::Point;  // inherit constructors

  Point& rotate(F angle, Point center) {
    Vector<2, F> v = (*this) - center;
    v.rotate(angle);
    Point p = center + v;
    (*this)[0] = p[0];
    (*this)[1] = p[1];

    return *this;
  }

  Point& rotate(F angle) {
    Vector<2, F> v;
    v[0] = (*this)[0];
    v[1] = (*this)[1];
    v.rotate(angle);
    (*this)[0] = v[0];
    (*this)[1] = v[1];

    return *this;
  }
};
}  // Geometry

#ifdef TEST_CASE
namespace Catch {

template <int N, typename F> struct StringMaker<Geometry::Vector<N, F>> {
  static std::string convert(const Geometry::Vector<N, F>& p) {
    std::ostringstream oss;
    oss << "<" << p[0];
    for (int i = 1; i < N; ++i)
      oss << ", " << p[i];
    oss << ">";
    return oss.str();
  }
};

template <int N, typename F> struct StringMaker<Geometry::Point<N, F>> {
  static std::string convert(const Geometry::Point<N, F>& p) {
    std::ostringstream oss;
    oss << "(" << p[0];
    for (int i = 1; i < N; ++i)
      oss << ", " << p[i];
    oss << ")";
    return oss.str();
  }
};
}
#endif
