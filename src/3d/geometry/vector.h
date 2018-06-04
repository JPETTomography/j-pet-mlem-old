#pragma once

#if !__CUDACC__
#include "util/json.h"
#include <istream>
#include "util/read.h"
#include <ostream>
#endif

#include "util/cuda/compat.h"
#include "2d/geometry/vector.h"

namespace PET3D {

/// 3D vector
template <typename FType> struct Vector {

  static Vector from_euler_angles(FType phi, FType theta) {
    FType r_xy = std::sin(theta);
    return Vector(r_xy * std::cos(phi), r_xy * std::sin(phi), std::cos(theta));
  }

  using F = FType;

  _ Vector(F x, F y, F z) : x(x), y(y), z(z) {}
  _ Vector() = default;

  F x, y, z;

#if !__CUDACC__
  /// construct Vector from json
  Vector(const json& j) : x(j[0]), y(j[1]), z(j[2]) {}
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

  _ Vector& operator*=(const Vector& v) {
    x *= v.x;
    y *= v.y;
    z *= v.z;
    return *this;
  }

  Vector& operator/=(FType s) {
    x /= s;
    y /= s;
    z /= s;
    return *this;
  }

  _ bool operator!=(const Vector& v) const {
    return x != v.x || y != v.y || z != v.z;
  }

  _ bool operator==(const Vector& v) const {
    return x == v.x && y == v.y && z == v.z;
  }

  _ F length2() const { return x * x + y * y + z * z; }

  _ F length() const { return compat::sqrt(length2()); }

  _ PET2D::Vector<FType> xy() const { return PET2D::Vector<FType>(x, y); }

  _ Vector normalized() const {
    Vector res(*this);
    res /= length();
    return res;
  }

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

  _ Vector operator-() const { return Vector(-x, -y, -z); }

  _ Vector operator*(FType rhs) {
    Vector vec(*this);
    vec *= rhs;
    return vec;
  }

  _ friend Vector operator*(FType lhs, const Vector& rhs) {
    Vector vec(rhs);
    vec *= lhs;
    return vec;
  }

  _ FType dot(const Vector& rhs) const {
    return x * rhs.x + y * rhs.y + z * rhs.z;
  }

  _ Vector cross(const Vector& rhs) const {
    return Vector(-z * rhs.y + y * rhs.z,  //
                  z * rhs.x - x * rhs.z,   //
                  -y * rhs.x + x * rhs.y);
  }

#if !__CUDACC__
  // serialize point to json
  operator json() const { return json{ x, y, z }; }

  friend std::ostream& operator<<(std::ostream& out, const Vector& vec) {
    out << vec.x << ' ' << vec.y << ' ' << vec.z;
    return out;
  }
#endif

  static Vector e_x() { return Vector(1, 0, 0); }
  static Vector e_y() { return Vector(0, 1, 0); }
  static Vector e_z() { return Vector(0, 0, 1); }
};
}

#ifdef TEST_CASE
namespace Catch {
template <typename FType> struct StringMaker<PET3D::Vector<FType>> {
  static std::string convert(const PET3D::Vector<FType>& vec) {
    std::ostringstream oss;
    oss << "<" << vec.x << ", " << vec.y << ", " << vec.z << ">";
    return oss.str();
  }
};
}
#endif
