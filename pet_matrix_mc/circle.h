#pragma once

#include "point.h"
#include "event.h"

// produces secant angles circle/line intersection as a equation system solution
// see /math/secant.nb
template <typename F = double>
class circle {
public:
  circle(F radious)
  : r(radious)
  , r2(radious*radious) {}

  typedef F angle_type;
  typedef point<F> point_type;
  typedef event<F> event_type;
  typedef std::pair<point_type, point_type> secant_type;

  secant_type secant(event_type &e) {
#if GENERIC_LINE_EQUATION
    auto sq  = sqrt( e.b2 * ( -(e.c*e.c) + e.a2_b2 * r2 ) );
    auto asq = e.a*sq;

    return std::make_pair(
      point_type( (e.ac - sq) / e.a2_b2, (e.b2c + asq) / e.b_a2_b2 ),
      point_type( (e.ac + sq) / e.a2_b2, (e.b2c - asq) / e.b_a2_b2 )
    );
#else
    auto sq  = sqrt(-e.b2 + e.m2*r2 + r2);
    auto msq = e.m*sq;

    if (e.flip) return std::make_pair(
      point_type( (e.b - msq) / e.m2p1, (-sq - e.bm ) / e.m2p1 ),
      point_type( (e.b + msq) / e.m2p1, ( sq - e.bm ) / e.m2p1 )
    );
    return std::make_pair(
      point_type( (-sq - e.bm ) / e.m2p1, (e.b - msq) / e.m2p1 ),
      point_type( ( sq - e.bm ) / e.m2p1, (e.b + msq) / e.m2p1 )
    );
#endif
  }

  std::pair<angle_type, angle_type>
  secant_angles(event_type &e) {
    auto s = secant(e);
    return std::make_pair( std::atan2(s.first.y,  s.first.x),
                           std::atan2(s.second.y, s.second.x) );
  }

  F radious()  const { return r;  }
  F radious2() const { return r2; }
private:
  const F r;
  const F r2;
};
