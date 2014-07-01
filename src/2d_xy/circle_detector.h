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

  CircleDetector(F radius)
      : Circle<F>(radius), center(), svg_class("detector") {}

  // this is for compatibility with square detector
  CircleDetector(F w, F h, F d __attribute__((unused)) = F())
      : Super(w / 2), center(), svg_class("detector") {
    if (w != h)
      throw("circle detector height and width must be equal");
  }

  static F default_height_for_width(const F w) { return w; }

  CircleDetector(F radius, Point center) : Super(radius), center(center) {}

  CircleDetector& rotate(Angle phi) {
    center.rotate(phi);
    return *this;
  }

  CircleDetector& operator+=(Point t) {
    center += t;
    return *this;
  }

  F max_distance() { return center.length() + this->radius(); }

  Point center;

  Intersections intersections(typename Super::Event e) {
    auto intersections = this->secant(e - center);
    for (auto& p : intersections) {
      p += center;
    }
    return intersections;
  }

  const char* svg_class;

  friend svg_ostream<F>& operator<<(svg_ostream<F>& svg, CircleDetector& cd) {
    svg << "<circle class=\"" << cd.svg_class << "\" cx=\"" << cd.center.x
        << "\" cy=\"" << cd.center.y << "\" r=\"" << cd.radius() << "\"/>"
        << std::endl;
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
