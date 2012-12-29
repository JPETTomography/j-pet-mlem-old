#pragma once

#include "point.h"
#include "event.h"
#include "svg_ostream.h"

// produces secant angles circle/line intersection as a equation system solution
// see /math/secant.nb
template <typename F = double>
class circle {
public:
  circle(F radius)
  : r(radius)
  , r2(radius*radius) {}

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


  F angle(point_type p) {
    return std::atan2(p.y, p.x);
  }

  std::pair<angle_type, angle_type>
  secant_angles(event_type &e) {
    auto s = secant(e);
    return std::make_pair( angle(s.first), angle(s.second));
  }


  int section(F angle,int n_detectors) {
    // converting angles to [0,2 Pi) interval
    F  normalised_angle=angle>0?angle:(F)2.0*M_PI+angle;
    return
      static_cast<int>(round( normalised_angle *n_detectors*INV_TWO_PI ))
                       % n_detectors;
  }

  std::pair<int, int>
  secant_sections(event_type &e, size_t n_detectors) {
    auto sa = secant_angles(e);

    return std::make_pair(section(sa.first,n_detectors),
                          section(sa.second,n_detectors));
  }

  F radius()  const { return r;  }
  F radius2() const { return r2; }

  friend svg_ostream<F> & operator << (svg_ostream<F> &svg, circle &c) {
    svg << "<circle cx=\"0\" cy=\"0\" r=\"" << c.r << "\"/>" << std::endl;
    return svg;
  }

private:
  const F r;
  const F r2;
  static const F TWO_PI;
  static const F INV_TWO_PI;
};

template<typename F> const F circle<F>::TWO_PI=(F)2.0*M_PI;
template<typename F> const F circle<F>::INV_TWO_PI=(F)1.0/TWO_PI;
