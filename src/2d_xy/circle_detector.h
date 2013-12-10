#pragma once

#include "geometry/point.h"
#include "geometry/circle.h"
#include "util/svg_ostream.h"
#include "util/array.h"

template <typename FType = double> class CircleDetector : Circle<FType> {
 public:
  typedef FType F;
  typedef Circle<F> Super;
  typedef F Angle;
  typedef ::Point<F> Point;
  typedef Array<2, Point> Intersections;
  typedef typename Super::Event Event;

  CircleDetector(F radius) : Circle<F>(radius), center() {}

  // this is for compatibility with square detector
  CircleDetector(F diameter, F height __attribute__((unused)))
      : Super(diameter / 2.), center() {}

  CircleDetector(F radius, Point center) : Super(radius), center(center) {}

  CircleDetector& rotate(Angle phi) {
    center.rotate(phi);
    return *this;
  }

  CircleDetector rotated(Angle phi) {
    return CircleDetector(this->radius(), center.rotated(phi));
  }

  CircleDetector& operator+=(Point t) {
    center += t;
    return *this;
  }

  F max_distance() { return center.length() + this->radius(); }

  Point center;

  Intersections intersections(typename Super::Event e) {
    return this->secant(e - center);
  }

  friend svg_ostream<F>& operator<<(svg_ostream<F>& svg, CircleDetector& cd) {
    svg << "<circle class=\"detector\" cx=\"" << cd.center.x << "\" cy=\""
        << cd.center.y << "\" r=\"" << cd.radius() << "\"/>" << std::endl;
    return svg;
  }

  friend std::ostream& operator<<(std::ostream& out, CircleDetector& cd) {
    out << "circle (" << cd.center.x << ", " << cd.center.y << ") radius "
        << cd.radius() << std::endl;
    out << std::flush;
    return out;
  }

 private:
  CircleDetector() {}
};
