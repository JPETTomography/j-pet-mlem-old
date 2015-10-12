#pragma once

#if !__CUDACC__
#include "util/json.h"
#include <istream>
#include "util/read.h"
#include <ostream>
#endif

#include "util/cuda/compat.h"

namespace PET2D {

/// 2D vector
template <typename FType> struct Vector {
  using F = FType;

  _ Vector(F x, F y) : x(x), y(y) {}
  _ Vector() = default;

  F x, y;

#if !__CUDACC__
  /// Construct Vector from JSON.
  Vector(const json& j) : x(j[0]), y(j[1]) {}

  /// Constructs Vector from stream.
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

  _ Vector& operator/=(F s) {
    x /= s;
    y /= s;
    return *this;
  }

  _ Vector& normalize() {
    F length = this->length();
    (*this /= length);
    return *this;
  }

  _ Vector normalized() {
    Vector res(*this);
    res.normalize();
    return res;
  }
  _ bool operator!=(const Vector& p) const { return x != p.x || y != p.y; }

  _ bool operator==(const Vector& p) const { return x == p.x && y == p.y; }

  _ F length2() const { return x * x + y * y; }

  _ F length() const { return compat::sqrt(x * x + y * y); }

  /// Rotate Vector around \f$ (0, 0) \f$ with given angle.
  ////
  /// \note
  /// I know it is bad idea to count all over again
  /// \c sin/cos for given Vector, but this will be used
  /// only for initialization.
  _ Vector& rotate(F phi) {
    F sin_phi = compat::sin(phi);
    F cos_phi = compat::cos(phi);
    F tx = x * cos_phi - y * sin_phi;
    F ty = x * sin_phi + y * cos_phi;
    x = tx;
    y = ty;
    return *this;
  }

  _ Vector rotated(F phi) const {
    Vector tmp(*this);
    return tmp.rotate(phi);
  }

  /// Return counter-clockwise perpendicular vector.
  _ Vector perpendicular() const { return Vector(-y, x); }

  /// Return clockwise perpendicular vector.
  _ Vector cw_perpendicular() const { return Vector(y, -x); }

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
  // serialize point to json
  operator json() const { return json{ x, y }; }

  friend std::ostream& operator<<(std::ostream& out, const Vector& vec) {
    out << vec.x << " " << vec.y;
    return out;
  }
#endif
};

}  // PET2D

#ifdef TEST_CASE
namespace Catch {
template <typename FType> struct StringMaker<PET2D::Vector<FType>> {
  static std::string convert(const PET2D::Vector<FType>& vec) {
    std::ostringstream oss;
    oss << "<" << vec.x << ", " << vec.y << ">";
    return oss.str();
  }
};
}
#endif
