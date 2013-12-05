#pragma once
#include <iostream>

#include "point.h"
#include "2d_xy/event.h"
#include "util/svg_ostream.h"

template <typename FType = double>
class Polygon : public std::vector<Point<FType>> {
 public:
  typedef FType F;
  typedef ::Point<F> Point;
  typedef ::Event<F> Event;
  typedef std::initializer_list<Point> Intersections;
  typedef ::svg_ostream<F> svg_ostream;

  // tests for intersection with generic form line equation
  bool intersects(Event& e) {
    auto p1 = this->back();
    auto v1 = e(p1);
    for (auto it = this->begin(); it != this->end(); ++it) {
      auto p2 = *it;
      auto v2 = e(p2);
      if (v1 * v2 <= 0.)
        return true;
      v1 = v2;
    }
    return false;
  }

  Intersections intersections(Event& e) {
    auto p1 = this->back();
    auto v1 = e(p1);
    Point first, second;
    auto count = 0;
    for (auto it = this->begin(); it != this->end(); ++it) {
      auto p2 = *it;
      auto v2 = e(p2);
      if (v2 == 0.) {
        // v2 is crossing point
        (!(count++) ? first : second) = p2;
      } else if (v1 * v2 < 0.) {
        // calculate intersection
        auto m = e.a * (p1.x - p2.x) + e.b * (p1.y - p2.y);
        (!(count++) ? first : second) = Point(
            (e.c * (p1.x - p2.x) + e.b * (p2.x * p1.y - p1.x * p2.y)) / m,
            (e.c * (p1.y - p2.y) + e.a * (p1.x * p2.y - p2.x * p1.y)) / m);
      }
      if (count == 2)

        return { first, second };
      v1 = v2;
      p1 = p2;
    }
    if (count)
      return { first };
    return {};
  }

  friend svg_ostream& operator<<(svg_ostream& svg, Polygon& pg) {
    svg << "<polygon points=\"";
    for (auto it = pg.begin(); it != pg.end(); ++it) {
      auto p = *it;
      svg << p.x << ' ' << p.y << ' ';
    }
    svg << "\"/>" << std::endl;
    return svg;
  }

  friend std::ostream& operator<<(std::ostream& out, Polygon& pg) {
    out << "polygon points \n";
    for (auto it = pg.begin(); it != pg.end(); ++it) {
      auto p = *it;
      out << p.x << ' ' << p.y << "\n";
    }
    out << std::flush;
    return out;
  }
};
