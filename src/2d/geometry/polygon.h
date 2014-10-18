#pragma once

#include <iostream>
#include <cmath>

#include "point.h"
#include "2d/barrel/event.h"
#include "util/svg_ostream.h"
#include "util/array.h"

namespace PET2D {

/// Closed polygon described by set of points
template <std::size_t NumPoints, typename FType = double>
class Polygon : public Array<NumPoints, Point<FType>> {
 public:
  typedef FType F;
  typedef F Angle;
  typedef PET2D::Point<F> Point;
  typedef Barrel::Event<F> Event;
  typedef Array<2, Point> Intersections;
  typedef ::svg_ostream<F> svg_ostream;

  Polygon& rotate(Angle phi) {
    for (auto& p : *this) {
      p.rotate(phi);
    }
    return *this;
  }

  Polygon& operator+=(Point t) {
    for (auto& p : *this) {
      p += t;
    }
    return *this;
  }

  // tests for intersection with generic form line equation
  bool intersects(Event& e) {
    auto p1 = this->back();
    auto v1 = e(p1);
    for (auto& p2 : *this) {
      auto v2 = e(p2);
      if (v1 * v2 <= 0)
        return true;
      v1 = v2;
    }
    return false;
  }

  F max_distance() {
    F distance = 0;
    for (auto& p : *this) {
      distance = std::max(distance, p.length());
    }
    return distance;
  }

  Intersections intersections(Event& e) {
    auto p1 = this->back();
    auto v1 = e(p1);
    Intersections intersections;
    for (auto& p2 : *this) {
      auto v2 = e(p2);
      if (v2 == 0) {
        // v2 is crossing point
        intersections.push_back(p2);
      } else if (v1 * v2 < 0) {
        // calculate intersection
        auto m = e.a * (p1.x - p2.x) + e.b * (p1.y - p2.y);
        intersections.push_back(Point(
            (e.c * (p1.x - p2.x) + e.b * (p2.x * p1.y - p1.x * p2.y)) / m,
            (e.c * (p1.y - p2.y) + e.a * (p1.x * p2.y - p2.x * p1.y)) / m));
      }
      if (intersections.full())
        return intersections;
      v1 = v2;
      p1 = p2;
    }
    return intersections;
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
}  // PET2D
