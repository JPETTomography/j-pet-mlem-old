#pragma once

#include "point.h"

#define GENERIC_LINE_EQUATION 1

template<typename F = double>
struct event {
  typedef point<F> point_type;

  event(F a_x, F a_y, F a_phi)
  : x(a_x)
  , y(a_y)
  , phi(a_phi)
  {
#if GENERIC_LINE_EQUATION
    // get line equation coefficients
    // a x + b y == c
    a =  std::sin(phi);
    b = -std::cos(phi);
    c = a*x + b*y;

    // helper variables
    b2      = b*b;
    b2c     = b2*c;
    ac      = a*c;
    a2_b2   = a*a + b2;
    b_a2_b2 = b * a2_b2;
#else
    // using normal line equation
    // switch coordinates when close to M_PI_2 to get accurate results
    flip = phi > M_PI_4 && phi < 3.0 * M_PI_4;
    m    = !flip ? tan(phi) : tan(phi - M_PI_2);
    b    = !flip ? (y - m*x) : (x - m*y);
    m2   = m*m;
    m2p1 = m2 + 1.;
    b2   = m*m;
    bm   = b*m;
#endif
  }

  event(point_type p, F phi)
  : event(p.x, p.y, phi) {}

  F operator () (const point_type &p) {
#if GENERIC_LINE_EQUATION
    return a*p.x + b*p.y - c;
#else
    return !flip ? (m*p.x + b - p.y) : (m*p.y + b - p.x);
#endif
  }

  F x, y;
  F phi;

  // line equation coefficients
#if GENERIC_LINE_EQUATION
  F a, b, c;
  // precalculated variables
  F b2, b2c, ac, a2_b2, b_a2_b2;
#else
  F m, b;
  // precalculated variables
  bool flip;
  F m2, m2p1, b2, bm;
#endif
};
