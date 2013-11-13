#pragma once

#include "geometry/point.h"

template <typename FType = double> class CircleDetector {
 public:
  typedef FType F;
  typedef F Angle;
  typedef ::Point<FType> Point;

  CircleDetector(F radius) : radius(radius), center() {}

  CircleDetector(F radius, Point center) : radius(radius), center(center) {}

  CircleDetector& rotate(Angle phi __attribute__((unused))) { return *this; }

  CircleDetector rotated(Angle phi __attribute__((unused))) {
    SquareDetector r(radius, center);
    return r;
  }

  CircleDetector& operator+=(Point t) {
    center += t;
    return *this;
  }

  F radius;
  Point center;

 private:
  CircleDetector() {}
};
