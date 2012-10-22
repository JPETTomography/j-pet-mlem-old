#pragma once

#include "point.h"

#define SECANT_GENERIC_LINE_EQUATION 0

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
  typedef std::pair<point_type, point_type> secant_type;

  secant_type secant(point_type &p, angle_type phi) {
#if SECANT_GENERIC_LINE_EQUATION
    // using generic line equation
    auto sin_phi = std::sin(phi);
    auto cos_phi = std::cos(phi);

    // get line equation coefficients
    // a x + b y == c
    auto a =  sin_phi;
    auto b = -cos_phi;
    auto c = p.x * sin_phi - p.y * cos_phi;

    // intersection points
    auto b2      = b*b;
    auto b2c     = b2*c;
    auto ac      = a*c;
    auto a2_b2   = a*a + b2;
    auto b_a2_b2 = b * a2_b2;
    auto sq      = sqrt( b2 * ( -(c*c) + a2_b2 * r2 ) );
    auto asq     = a*sq;

    return std::make_pair(
      point_type( (ac - sq) / a2_b2, (b2c + asq) / b_a2_b2 ),
      point_type( (ac + sq) / a2_b2, (b2c - asq) / b_a2_b2 )
    );
#else
    // using normal line equation
    // switch coordinates when close to M_PI_2 to get accurate results
    bool degen = phi > M_PI_4 && phi < 3.0 * M_PI_4;
    auto m     = !degen ? tan(phi) : tan(phi - M_PI_2);
    auto b     = !degen ? (p.y - m*p.x) : (p.x - m*p.y);
    auto m2    = m*m;
    auto m2p1  = m2 + 1.;
    auto b2    = m*m;
    auto sq    = sqrt(-b2 + m2*r2 + r2);

    // switch coordinates on degen
    if (degen) return std::make_pair(
      point_type( (b - m*sq) / m2p1, (-sq - b*m ) / m2p1 ),
      point_type( (b + m*sq) / m2p1, ( sq - b*m ) / m2p1 )
    );
    return std::make_pair(
      point_type( (-sq - b*m ) / m2p1, (b - m*sq) / m2p1 ),
      point_type( ( sq - b*m ) / m2p1, (b + m*sq) / m2p1 )
    );
#endif
  }

  std::pair<angle_type, angle_type>
  secant_angles(point_type &p, angle_type phi) {
    auto s = secant(p, phi);
    return std::make_pair( std::atan2(s.first.y,  s.first.x),
                           std::atan2(s.second.y, s.second.x) );
  }

  F radious()  const { return r;  }
  F radious2() const { return r2; }
private:
  const F r;
  const F r2;
};
