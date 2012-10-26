#pragma once

#include "point.h"
#include "event.h"
#include "svg_ostream.h"

template <typename F = double>
class polygon : public std::vector<point<F>> {
public:
  typedef point<F> point_type;
  typedef event<F> event_type;
  typedef std::vector<point_type> intersections_type;

  // tests for intersection with generic form line equation
  bool intersects(event_type &e) {
    auto p1 = this->back();
    auto v1 = e(p1);
    for(auto p2: *this) {
      auto v2 = e(p2);
      if (v1 * v2 <= 0.) return true;
      v1 = v2;
    }
    return false;
  }

  intersections_type
  intersections(event_type &e) {
    auto first = false;
    auto p1 = this->back();
    auto v1 = e(p1);
    intersections_type r;
    for(auto p2: *this) {
      auto v2 = e(p2);
      if (v2 == 0.) {
        // v2 is crossing point
        r.push_back(p2);
        if (r.size() == 2) return r;
      } else if (v1 * v2 < 0.) {
        // calculate intersection
#if GENERIC_LINE_EQUATION
        auto m = e.a*(p1.x - p2.x) + e.b*(p1.y - p2.y);
        r.push_back({
          ( e.c*(p1.x - p2.x) + e.b*(p2.x*p1.y - p1.x*p2.y) ) / m,
          ( e.c*(p1.y - p2.y) + e.a*(p1.x*p2.y - p2.x*p1.y) ) / m
        });
#endif
        if (r.size() == 2) return r;
      }
      v1 = v2;
      p1 = p2;
    }
    return r;
  }

  friend svg_ostream<F> & operator << (svg_ostream<F> &svg, polygon &pg) {
    svg << "<polygon points=\"";
    for(auto p: pg) {
      svg << p.x << ' ' << p.y << ' ';
    }
    svg << "\"/>" << std::endl;
    return svg;
  }
};
