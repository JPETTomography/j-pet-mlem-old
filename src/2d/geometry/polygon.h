#pragma once

#include "point.h"
#include "2d/barrel/event.h"
#include "util/array.h"
#include "util/cuda/compat.h"
#if !__CUDACC__
#include "util/svg_ostream.h"
#endif

namespace PET2D {

/// Closed polygon described by set of points
template <std::size_t NumPoints, typename FType = double>
class Polygon : public util::array<NumPoints, Point<FType>> {
 public:
  using F = FType;
  using Angle = F;
  using Point = PET2D::Point<F>;
  using Event = Barrel::Event<F>;
  using Intersections = util::array<2, Point>;

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

  /// Returns center point of the polygon
  Point center() const {
    Point center_point(0, 0);
    for (auto& p : *this) {
      center_point += p;
    }
    center_point.x /= this->size();
    center_point.y /= this->size();
    return center_point;
  }

  // tests for intersection with generic form line equation
  bool intersects(Event& e) const {
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
      distance = compat::max(distance, p.length());
    }
    return distance;
  }

  _ Intersections intersections(Event& e) const {
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
        intersections.emplace_back(
            (e.c * (p1.x - p2.x) + e.b * (p2.x * p1.y - p1.x * p2.y)) / m,
            (e.c * (p1.y - p2.y) + e.a * (p1.x * p2.y - p2.x * p1.y)) / m);
      }
      if (intersections.full())
        return intersections;
      v1 = v2;
      p1 = p2;
    }
    return intersections;
  }

#if !__CUDACC__
  using svg_ostream = util::svg_ostream<F>;

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
#endif
};
}  // PET2D
