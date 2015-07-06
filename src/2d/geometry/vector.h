#pragma once

#if !__CUDACC__
#include <ostream>
#endif

#include "util/cuda/compat.h"
#include "util/read.h"

namespace PET2D {

/// 2D Vector with given coordinates
template <typename FType> struct Vector {
  using F = FType;

  _ Vector(F x, F y) : x(x), y(y) {}
  _ Vector() = default;

  F x, y;

#if !__CUDACC__
  /// constructs Vector from stream
  Vector(std::istream& in) : x(util::read<F>(in)), y(util::read<F>(in)) {}
#endif

  _ Vector& operator+=(const Vector& p) {
    x += p.x;
    y += p.y;
    return *this;
  }

  _ Vector& operator-=(const Vector& p) {
    x -= p.x;
    y -= p.y;
    return *this;
  }

  _ Vector& operator*=(F s) {
    x *= s;
    y *= s;
    return *this;
  }

  Vector& operator/=(F s) {
    x /= s;
    y /= s;
    return *this;
  }

  Vector& normalize() {
    F length = this->length();
    (*this /= length);
    return *this;
  }

  Vector normalized() {
    Vector res(*this);
    res.normalize();
    return res;
  }
  _ bool operator!=(const Vector& p) const { return x != p.x || y != p.y; }

  _ bool operator==(const Vector& p) const { return x == p.x && y == p.y; }

  _ F length2() const { return x * x + y * y; }

  _ F length() const { return compat::sqrt(x * x + y * y); }

  /// Rotate Vector around (0, 0) with given angle

  /// \note
  /// I know it is bad idea to count all over again
  /// \c sin/cos for given Vector, but this will be used
  /// only for initialization.
  Vector& rotate(F phi) {
    F sin_phi = compat::sin(phi);
    F cos_phi = compat::cos(phi);
    F tx = x * cos_phi - y * sin_phi;
    F ty = x * sin_phi + y * cos_phi;
    x = tx;
    y = ty;
    return *this;
  }

  Vector rotated(F phi) const {
    Vector tmp(*this);
    return tmp.rotate(phi);
  }

  Vector perpendicular() const { return Vector(-y, x); }

  _ Vector operator+(const Vector& rhs) const {
    Vector vec(*this);
    vec += rhs;
    return vec;
  }

  _ Vector operator-(const Vector& rhs) const {
    Vector vec(*this);
    vec -= rhs;
    return vec;
  }

  _ F dot(const Vector& rhs) const { return x * rhs.x + y * rhs.y; }

#if !__CUDACC__
  friend std::ostream& operator<<(std::ostream& out, const Vector& vec) {
    out << vec.x << " " << vec.y;
    return out;
  }
#endif
};
}
