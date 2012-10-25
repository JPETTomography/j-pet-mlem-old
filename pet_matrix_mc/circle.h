#pragma once

#include "point.h"
#include "event.h"
#include "svg_ostream.h"

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
#if DEBUG
    std::cout
      << " x=" << e.x
      << " y=" << e.y
      << " phi=" << e.phi
      << std::endl;
#endif
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
#if DEBUG
    std::cout
      << " x1=" << s.first.x
      << " y1=" << s.first.y
      << " x2=" << s.second.x
      << " y2=" << s.second.y
      << std::endl;
#endif
    return std::make_pair( std::atan2(s.first.y,  s.first.x),
                           std::atan2(s.second.y, s.second.x) );
  }

  std::pair<size_t, size_t>
  secant_sections(event_type &e, size_t n_detectors) {
    auto sa = secant_angles(e);
#if DEBUG
    std::cout
      << "a1=" << sa.first
      << "a2=" << sa.second
      << std::endl;
#endif
    return std::make_pair(
      static_cast<int>( round( sa.first  * n_detectors / (2. * M_PI) ) ) % n_detectors,
      static_cast<int>( round( sa.second * n_detectors / (2. * M_PI) ) ) % n_detectors
    );
  }

  F radious()  const { return r;  }
  F radious2() const { return r2; }

  friend svg_ostream<F> & operator << (svg_ostream<F> &svg, circle &c) {
    svg << "<circle cx=\"0\" cy=\"0\" r=\"" << c.r << "\" stroke=\"green\" fill=\"transparent\" stroke-width=\".01\"/>" << std::endl;
    return svg;
  }

private:
  const F r;
  const F r2;
};
