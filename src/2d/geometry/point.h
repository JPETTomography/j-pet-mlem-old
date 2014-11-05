#pragma once

#include <cmath>

#include "util/cuda/compat.h"

#include "pixel.h"

namespace PET2D {

/// 2D point with given coordinates
template <typename FType = double, typename SType = int> struct Point {
  using F = FType;
  using S = SType;
  using Pixel = PET2D::Pixel<S>;

  _ Point(F x, F y) : x(x), y(y) {}
  _ Point() = default;

  F x, y;

#if !__CUDACC__
  /// constructs Point from stream
  Point(std::istream& in) : x(util::read<F>(in)), y(util::read<F>(in)) {}
#endif

  _ Point operator+(const Point& p) const { return Point(x + p.x, y + p.y); }

  _ Point operator-(const Point& p) const { return Point(x - p.x, y - p.y); }

  _ Point& operator+=(const Point& p) {
    x += p.x;
    y += p.y;
    return *this;
  }

  _ Point& operator-=(const Point& p) {
    x -= p.x;
    y -= p.y;
    return *this;
  }

  _ bool operator!=(const Point& p) const { return x != p.x || y != p.y; }

  _ bool operator==(const Point& p) const { return x == p.x && y == p.y; }

  _ bool operator<(const Point& p) const {
    return y < p.y || (y == p.y && x < p.x);
  }

  F length2() const { return x * x + y * y; }

  F length() const { return compat::sqrt(x * x + y * y); }

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

  F nearest_distance(const Point& p1, const Point& p2) const {
    return compat::min((p1 - *this).length(), (p2 - *this).length());
  }

  Pixel pixel(F pixel_size, S pixel_count_2) {
    return Pixel(static_cast<S>(std::floor(x / pixel_size + pixel_count_2)),
                 static_cast<S>(std::floor(y / pixel_size + pixel_count_2)));
  }
};

/// Single point source
template <typename FType = double, typename SType = int>
struct PointSource : public Point<FType, SType> {
  using F = FType;
  using S = SType;
  using Point = PET2D::Point<F, S>;

  const F intensity;

  PointSource(Point p, F intensity) : Point::Point(p), intensity(intensity) {}

#if !__CUDACC__
  /// constructs point source from stream
  PointSource(std::istream& in)
      : Point::Point(in), intensity(util::read<F>(in)) {}
#endif
};
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
