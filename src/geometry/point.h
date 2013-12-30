#pragma once

#include <cmath>

template <typename FType = double> struct Point {
  typedef FType F;

  Point() : x(static_cast<F>(0)), y(static_cast<F>(0)) {}
  Point(F x, F y) : x(x), y(y) {}
  Point(std::pair<F, F> p) : x(p.first), y(p.second) {}

  F x, y;

  // I know it is bad idea to count all over again
  // sin/cos for given point, but this will be used
  // only for initialization.

  Point& rotate(F phi) {
    auto sin_phi = std::sin(phi);
    auto cos_phi = std::cos(phi);
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

  Point operator+(const Point& p) const {
    return { x + p.x, y + p.y };
  }

  Point operator-(const Point& p) const {
    return { x - p.x, y - p.y };
  }

  Point& operator+=(const Point& p) {
    x += p.x;
    y += p.y;
    return *this;
  }

  F length() const { return std::sqrt(x * x + y * y); }

  F nearest_distance(const Point& p1, const Point& p2) const {
    return std::min((p1 - *this).length(), (p2 - *this).length());
  }
};

template <typename F> F deg(F rad) { return rad * 180 / M_PI; }
template <typename F> F rad(F deg) { return deg * M_PI / 180; }