#pragma once

#if !__CUDACC__
#include "util/json.h"
#include <istream>
#include "util/read.h"
#include <ostream>
#endif

#include "util/cuda/compat.h"

#include "pixel.h"
#include "2d/geometry/vector.h"

namespace PET2D {

/// 2D point with given coordinates
template <typename FType> struct Point {
  using F = FType;
  using Vector = PET2D::Vector<FType>;

  _ Point(F x, F y) : x(x), y(y) {}
  _ Point() = default;

  F x, y;

#if !__CUDACC__
  /// construct Point from json
  Point(const json& j) : x(j[0]), y(j[1]) {}

  /// construct Point from stream
  Point(std::istream& in) : x(util::read<F>(in)), y(util::read<F>(in)) {}
#endif

  _ Point& operator+=(const Vector& v) {
    x += v.x;
    y += v.y;
    return *this;
  }

  _ Point& operator-=(const Vector& v) {
    x -= v.x;
    y -= v.y;
    return *this;
  }

  _ bool operator==(const Point& p) const { return x == p.x && y == p.y; }

  _ bool operator!=(const Point& p) const { return !this->operator==(p); }

  _ F distance_from_origin2() const { return x * x + y * y; }

  _ F distance_from_origin() const { return compat::sqrt(x * x + y * y); }

  /// Rotate point around (0, 0) with given angle

  /// \note
  /// I know it is bad idea to count all over again
  /// \c sin/cos for given point, but this will be used
  /// only for initialization.
  Point& rotate(F phi) {
    F sin_phi = compat::sin(phi);
    F cos_phi = compat::cos(phi);
    F tx = x * cos_phi - y * sin_phi;
    F ty = x * sin_phi + y * cos_phi;
    x = tx;
    y = ty;
    return *this;
  }

  Point rotated(F phi) const {
    Point tmp(*this);
    return tmp.rotate(phi);
  }

  _ F nearest_distance(const Point& p1, const Point& p2) const {
    return compat::min((p1 - *this).length(), (p2 - *this).length());
  }

  template <typename PType> _ PType pixel(F pixel_size, int pixel_count_2) {
    return PType(static_cast<typename PType::S>(
                     compat::floor(x / pixel_size + pixel_count_2)),
                 static_cast<typename PType::S>(
                     compat::floor(y / pixel_size + pixel_count_2)));
  }

  _ Vector as_vector() const { return Vector(x, y); }

  _ Point operator+(const Vector& rhs) const {
    Point p(*this);
    p += rhs;
    return p;
  }

  _ Point operator-(const Vector& rhs) const {
    Point p(*this);
    p -= rhs;
    return p;
  }

  _ Vector operator-(const Point& rhs) const {
    return Vector(x - rhs.x, y - rhs.y);
  }

  _ Point interpolate(const Point& end, F t) const {
    return Point(x * (1 - t) + end.x * t, y * (1 - t) + end.y * t);
  }

#if !__CUDACC__
  // serialize point to json
  operator json() const { return json{ x, y }; }

  friend std::ostream& operator<<(std::ostream& out, const Point& vec) {
    out << vec.x << " " << vec.y;
    return out;
  }
#endif
};

}  // PET2D

template <typename FType> FType deg(FType rad) {
  return rad * 180 / FType(M_PI);
}
template <typename FType> FType rad(FType deg) {
  return deg * FType(M_PI) / 180;
}

#ifdef TEST_CASE
namespace Catch {
template <typename FType> struct StringMaker<PET2D::Point<FType>> {
  static std::string convert(const PET2D::Point<FType>& p) {
    std::ostringstream oss;
    oss << "(" << p.x << ", " << p.y << ")";
    return oss.str();
  }
};
}
#endif
