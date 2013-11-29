#pragma once

#include "geometry/point.h"
#include "geometry/circle.h"

template <typename FType = double> class CircleDetector : Circle<FType> {
 public:
  typedef FType F;
  typedef F Angle;
  typedef ::Point<FType> Point;

  CircleDetector(F radius) : Circle<F>(radius), center() {}

  // this is for compatibility with square detector
  CircleDetector(F diameter, F diameter2) : Circle<F>(diameter / 2.), center() {
    if (diameter != diameter2)
      throw("circle detector width and height must be equal");
  }

  CircleDetector(F radius, Point center) : Circle<F>(radius), center(center) {}

  CircleDetector& rotate(Angle phi __attribute__((unused))) { return *this; }

  CircleDetector rotated(Angle phi __attribute__((unused))) {
    return CircleDetector(this->radius(), center);
  }

  CircleDetector& operator+=(Point t) {
    center += t;
    return *this;
  }

  Point center;

 private:
  CircleDetector() {}
};
