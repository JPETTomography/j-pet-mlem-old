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
    auto sq  = sqrt( e.b2 * ( -(e.c*e.c) + e.a2_b2 * r2 ) );
    auto asq = e.a*sq;

    return std::make_pair(
      point_type( (e.ac - sq) / e.a2_b2, (e.b2c + asq) / e.b_a2_b2 ),
      point_type( (e.ac + sq) / e.a2_b2, (e.b2c - asq) / e.b_a2_b2 )
    );
  }

  std::pair<angle_type, angle_type>
  secant_angles(event_type &e) {
    auto s = secant(e);
    return std::make_pair( std::atan2(s.first.y,  s.first.x),
                           std::atan2(s.second.y, s.second.x) );
  }

  std::pair<size_t, size_t>
  secant_sections(event_type &e, size_t n_detectors) {
    auto sa = secant_angles(e);
    auto angle_1=sa.first>0?sa.first:2*M_PI+sa.first;
    auto angle_2=sa.second>0?sa.second:2*M_PI+sa.second;

    
    return std::make_pair(
      static_cast<int>( floor( angle_1  * n_detectors / (2.0 * M_PI) ) ) % n_detectors,
      static_cast<int>( floor( angle_2 * n_detectors / (2.0 * M_PI) ) ) % n_detectors
    );
  }

  F radious()  const { return r;  }
  F radious2() const { return r2; }

  friend svg_ostream<F> & operator << (svg_ostream<F> &svg, circle &c) {
    svg << "<circle cx=\"0\" cy=\"0\" r=\"" << c.r << "\"/>" << std::endl;
    return svg;
  }

private:
  const F r;
  const F r2;
};
