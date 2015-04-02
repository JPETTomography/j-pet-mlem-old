#pragma once

#include "util/cuda/compat.h"
#include "util/read.h"

namespace PET3D {

/// 2D Vector with given coordinates
template <typename FType = double, typename SType = int> struct Vector {
  using F = FType;
  using S = SType;

  _ Vector(F x, F y, F z) : x(x), y(y), z(z) {}
  _ Vector() = default;

  F x, y, z;

#if !__CUDACC__
  /// constructs Vector from stream
  Vector(std::istream& in) : x(util::read<F>(in)), y(util::read<F>(in)) {}
#endif

  _ Vector& operator+=(const Vector& v) {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
  }

  _ Vector& operator-=(const Vector& v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
  }

  _ Vector& operator*=(FType s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
  }

  Vector& operator/=(FType s) {
    x /= s;
    y /= s;
    z /= s;
    return *this;
  }



  _ bool operator!=(const Vector& p) const { return x != p.x || y != p.y; }

  _ bool operator==(const Vector& p) const { return x == p.x && y == p.y; }

  _ bool operator<(const Vector& p) const {
    return y < p.y || (y == p.y && x < p.x);
  }

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
};

template <typename FType>
_ Vector<FType> operator+(const Vector<FType>& lhs, const Vector<FType>& rhs) {
  Vector<FType> vec(lhs);
  vec += rhs;
  return vec;
}

template <typename FType>
_ Vector<FType> operator-(const Vector<FType>& lhs, const Vector<FType>& rhs) {
  Vector<FType> vec(lhs);
  vec -= rhs;
  return vec;
}
}
