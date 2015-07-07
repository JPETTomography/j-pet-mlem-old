#pragma once

#if !__CUDACC__
#include <ostream>
#endif

#include "util/cuda/compat.h"

#include "3d/geometry/vector.h"
#include "2d/geometry/point.h"

namespace PET3D {

/// 3D point with given coordinates
template <typename FType> struct Point {
  using F = FType;
  using Vector = PET3D::Vector<FType>;
  using Point2D = PET2D::Point<F>;

  _ Point(F x, F y, F z) : x(x), y(y), z(z) {}
  _ Point() = default;

  F x, y, z;

#if !__CUDACC__
  /// constructs Point from stream
  Point(std::istream& in) : x(util::read<F>(in)), y(util::read<F>(in)) {}
#endif

  _ Point& operator+=(const Vector& v) {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
  }

  _ Point& operator-=(const Vector& v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
  }

  _ bool operator==(const Point& p) const {
    return x == p.x && y == p.y && z == p.z;
  }
  _ bool operator!=(const Point& p) const { return !this->operator==(p); }

  _ F distance_from_origin2() const { return as_vector().length2(); }

  _ F distance_from_origin() const { return as_vector().length(); }

  _ F nearest_distance(const Point& p1, const Point& p2) const {
    return compat::min((p1 - *this).length(), (p2 - *this).length());
  }

  _ Vector as_vector() const { return Vector(x, y, z); }

  _ static Point from_vector(const Vector& vec) {
    return Point(vec.x, vec.y, vec.z);
  }

  _ Point2D xy() const { return Point2D(x, y); }

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
    return Vector(x - rhs.x, y - rhs.y, z - rhs.z);
  }

  _ Point interpolate(const Point& end, F t) const {
    return Point(x * (1 - t) + end.x * t,
                 y * (1 - t) + end.y * t,
                 z * (1 - t) + end.z * t);
  }

#if !__CUDACC__
  friend std::ostream& operator<<(std::ostream& out, const Point& p) {
    out << p.x << " " << p.y << "  " << p.z;
    return out;
  }
#endif
};

/// Single point source
template <typename FType> struct PointSource : public Point<FType> {
  using F = FType;
  using Point = PET3D::Point<F>;

  const F intensity;

  PointSource(Point p, F intensity) : Point::Point(p), intensity(intensity) {}

#if !__CUDACC__
  /// constructs point source from stream
  PointSource(std::istream& in)
      : Point::Point(in), intensity(util::read<F>(in)) {}
#endif
};

}  // PET3D

#ifdef TEST_CASE
namespace Catch {
template <typename FType> struct StringMaker<PET3D::Point<FType>> {
  static std::string convert(const PET3D::Point<FType>& p) {
    std::ostringstream oss;
    oss << "(" p.x << " " << p.y << " " << p.z << ")";
    return oss.str();
  }
};
}
#endif
