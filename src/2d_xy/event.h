#pragma once

#include "geometry/point.h"

template <typename FType = double> struct Event {
  typedef FType F;
  typedef ::Point<F> Point;

  Event(F x_a, F y_a, F phi_a) : x(x_a), y(y_a), phi(phi_a) {
    // get line equation coefficients
    // a x + b y == c
    a = std::sin(phi);
    b = -std::cos(phi);
    c = a * x + b * y;

    // helper variables
    b2 = b * b;
    b2c = b2 * c;
    ac = a * c;
    a2_b2 = a * a + b2;
    b_a2_b2 = b * a2_b2;
  }

  Event(Point p, F phi) : Event(p.x, p.y, phi) {}

  // evaluates line equation side on given point
  // 0 means points lies on the line, -1 left, 1 right
  F operator()(const Point& p) { return a * p.x + b * p.y - c; }

  F x, y;
  F phi;

  // line equation coefficients
  F a, b, c;
  // precalculated variables
  F b2, b2c, ac, a2_b2, b_a2_b2;
};
