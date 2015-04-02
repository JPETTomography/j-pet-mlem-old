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



  _ bool operator!=(const Vector& v) const { return x != v.x || y != v.y || z != v.z; }

  _ bool operator==(const Vector& v) const { return x == v.x && y == v.y && z == v.z; }


  _ F length2() const { return x * x + y * y +z*z; }

  _ F length() const { return compat::sqrt(length2()); }


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
