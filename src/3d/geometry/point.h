#pragma once

#include "util/cuda/compat.h"

#include "3d/geometry/vector.h"

namespace PET3D {

/// 3D point with given coordinates
template <typename FType = double, typename SType = int> struct Point {
  using F = FType;
  using S = SType;
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
};

/// Single point source
template <typename FType = double, typename SType = int>
struct PointSource : public Point<FType, SType> {
  using F = FType;
  using S = SType;
  using Point = PET3D::Point<F, S>;

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

}  // PET2D



#ifdef TEST_CASE
namespace Catch {
template <typename FType> struct StringMaker<PET3D::Point<FType>> {
  static std::string convert(const PET3D::Point<FType>& p) {
    std::ostringstream oss;
    oss << "(" << p.x << ", " << p.y << p.z << ")";
    return oss.str();
  }
};
}
#endif
