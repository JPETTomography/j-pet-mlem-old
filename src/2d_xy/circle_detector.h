#pragma once

#include "geometry/point.h"
#include "geometry/circle.h"

template <typename FType = double> class CircleDetector : Circle<FType> {
 public:
  typedef FType F;
  typedef Circle<F> Super;
  typedef F Angle;
  typedef ::Point<F> Point;
  typedef std::initializer_list<Point> Intersections;
  typedef typename Super::Event Event;

  CircleDetector(F radius) : Circle<F>(radius), center() {}

  // this is for compatibility with square detector
  CircleDetector(F diameter, F diameter2) : Super(diameter / 2.), center() {
    if (diameter != diameter2)
      throw("circle detector width and height must be equal");
  }

  CircleDetector(F radius, Point center) : Super(radius), center(center) {}

  // rotation has no effect on circle (return itself)
  CircleDetector& rotate(Angle phi __attribute__((unused))) { return *this; }
  CircleDetector rotated(Angle phi __attribute__((unused))) { return *this; }

  CircleDetector& operator+=(Point t) {
    center += t;
    return *this;
  }

  Point center;

  Intersections intersections(typename Super::Event e) {
    return this->secant(e - center);
  }

 private:
  CircleDetector() {}
};
