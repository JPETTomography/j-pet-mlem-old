#pragma once

#include <ostream>

#include "util/cuda/compat.h"

#include "3d/geometry/vector.h"
#include "2d/geometry/point.h"

namespace PET3D {

/// 3D point with given coordinates
template <typename FType> struct Point {
  using F = FType;
  using Vector = PET3D::Vector<FType>;

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

  _ PET2D::Point<F> xy() const { return PET2D::Point<F>(x, y); }
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

template <typename F>
_ Point<F> operator+(const Point<F>& lhs, const Vector<F>& rhs) {
  Point<F> p(lhs);
  p += rhs;
  return p;
}

template <typename F>
_ Point<F> operator-(const Point<F>& lhs, const Vector<F>& rhs) {
  Point<F> p(lhs);
  p -= rhs;
  return p;
}

template <typename F>
_ Vector<F> operator-(const Point<F>& lhs, const Point<F>& rhs) {
  return Vector<F>(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
}

template <typename FType>
std::ostream& operator<<(std::ostream& out, const Point<FType>& p) {
  out << p.x << " " << p.y << "  " << p.z;
  return out;
}

template <typename F> _ Point<F> from_vector(const Vector<F>& vec) {
  return Point<F>(vec.x, vec.y, vec.z);
}

template <typename FType> _
Point<FType> interpolate(FType t, const Point<FType>& start, const Point<FType>& end) {
  return Point<FType>(start.x * (1 - t) + end.x * t,
                      start.y * (1 - t) + end.y * t,
                      start.z * (1 - t) + end.z * t);
}

}  // PET3D

#ifdef TEST_CASE
namespace Catch {
template <typename FType> struct StringMaker<PET3D::Point<FType>> {
  static std::string convert(const PET3D::Point<FType>& p) {
    std::ostringstream oss;
    oss << p.x << " " << p.y << " " << p.z;
    return oss.str();
  }
};
}

#endif
