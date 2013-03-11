#pragma once

#include <cmath>

template <typename FType = double> struct Point {
  typedef FType F;

  Point(F x_a  , F y_a ) : x(x_a), y(y_a) {}
  Point(std::pair<F, F> p) : x(p.first), y(p.second) {}

  F x, y;

  // I know it is bad idea to count all over again
  // sin/cos for given point, but this will be used
  // only for initialization.
  Point rotated(F phi) {
    auto sin_phi = std::sin(phi);
    auto cos_phi = std::cos(phi);
    return { x * cos_phi - y * sin_phi, x * sin_phi + y * cos_phi };
  }

  Point operator+(const Point& p) const {
    return { x + p.x, y + p.y };
  }

  Point operator-(const Point& p) const{
    return { x - p.x, y - p.y };
  }

  Point& operator+=(const Point & p) {
    x += p.x;
    y += p.y;
    return *this;
  }

  F length() { return std::sqrt(x * x + y * y); }
};

template <typename F> F deg(F rad) { return rad * 180 / M_PI; }
template <typename F> F rad(F deg) { return deg * M_PI / 180; }
