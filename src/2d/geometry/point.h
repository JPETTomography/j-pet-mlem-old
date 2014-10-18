#pragma once

#include <cmath>

#include "util/cuda/compat.h"

#include "pixel.h"

namespace PET2D {

template <typename FType = double, typename SType = int> struct Point {
  typedef FType F;
  typedef SType S;
  typedef PET2D::Pixel<S> Pixel;

  _ Point() : x(0), y(0) {}
  _ Point(F x, F y) : x(x), y(y) {}

  F x, y;

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

  // I know it is bad idea to count all over again
  // sin/cos for given point, but this will be used
  // only for initialization.
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
}  // PET2D

template <typename F> F deg(F rad) { return rad * 180 / F(M_PI); }
template <typename F> F rad(F deg) { return deg * F(M_PI) / 180; }

#ifdef TEST_CASE
namespace Catch {
template <typename FType> struct StringMaker</**/ ::Point<FType>> {
  static std::string convert(const ::Point<FType>& p) {
    std::ostringstream oss;
    oss << "(" << p.x << ", " << p.y << ")";
    return oss.str();
  }
};
}
#endif
