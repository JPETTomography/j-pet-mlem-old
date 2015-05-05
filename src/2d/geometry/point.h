#pragma once

#include "util/cuda/compat.h"
#include "pixel.h"

#include "2d/geometry/vector.h"

namespace PET2D {

/// 2D point with given coordinates
template <typename FType> struct Point {
  using F = FType;

  // using Pixel = PET2D::Pixel<S>;
  using Vector = PET2D::Vector<FType>;

  _ Point(F x, F y) : x(x), y(y) {}
  _ Point() = default;

  F x, y;

#if !__CUDACC__
  /// constructs Point from stream
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
};

/// Single point source
template <typename FType> struct PointSource : public Point<FType> {
  using F = FType;
  using Point = PET2D::Point<F>;

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
  return Vector<F>(lhs.x - rhs.x, lhs.y - rhs.y);
}

}  // PET2D

template <typename F> F deg(F rad) { return rad * 180 / F(M_PI); }
template <typename F> F rad(F deg) { return deg * F(M_PI) / 180; }

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
