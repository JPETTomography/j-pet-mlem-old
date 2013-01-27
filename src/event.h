#pragma once

#include "point.h"

template <typename F = double> struct event {
  typedef point<F> point_type;

  event(F a_x, F a_y, F a_phi) : x(a_x), y(a_y), phi(a_phi) {
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

  event(point_type p, F phi) : event(p.x, p.y, phi) {}

  // evaluates line equation side on given point
  // 0 means points lies on the line, -1 left, 1 right
  F operator()(const point_type& p) { return a * p.x + b * p.y - c; }

  F x, y;
  F phi;

  // line equation coefficients
  F a, b, c;
  // precalculated variables
  F b2, b2c, ac, a2_b2, b_a2_b2;
};
