#pragma once

#include "point.h"
#include "2d/barrel/event.h"
#include "util/array.h"
#include "util/cuda/compat.h"
#include "transformation.h"

#if !__CUDACC__
#include "util/json.h"
#include "util/svg_ostream.h"
#endif

namespace PET2D {

/// Closed polygon described by set of points
template <std::size_t NumPoints, typename FType>
class Polygon : public util::array<NumPoints, Point<FType>> {
 public:
  using F = FType;
  using Angle = F;
  using Point = PET2D::Point<F>;
  using Vector = PET2D::Vector<F>;
  using Event = Barrel::Event<F>;
  using Transformation = PET2D::Transformation<F>;
  using Intersections = util::array<2, Point>;

  Polygon& transform(Transformation tr) {
    for (auto& p : *this) {
      p = tr(p);
    }
    return *this;
  }

  Polygon& transformed(Transformation tr) {
    Polygon tmp(*this);
    return tmp.transform(tr);
  }

  Polygon& rotate(Angle phi) {
    for (auto& p : *this) {
      p.rotate(phi);
    }
    return *this;
  }

  Polygon& operator+=(Vector t) {
    for (auto& p : *this) {
      p += t;
    }
    return *this;
  }

  Polygon operator+(Vector t) const { return Polygon(*this) += t; }

  /// Returns center point of the polygon
  Point center() const {
    //  Point center_point(0, 0);
    FType c_x = 0.0;
    FType c_y = 0.0;
    for (auto& p : *this) {
      c_x += p.x;
      c_y += p.y;
    }
    c_x /= this->size();
    c_y /= this->size();
    return Point(c_x, c_y);
  }

  bool contains(const Point& p, F tol = 0) {
    int n = this->size();
    Point start = (*this)[n - 1];

    for (int i = 0; i < n; i++) {
      Point end = (*this)[i];
      Vector edge = end - start;
      if (edge.cw_perpendicular().dot(p - start) < -tol) {
        return false;
      }
      start = end;
    }
    return true;
  }

  // tests for intersection with generic form line equation
  bool intersects(Event& event) const {
    auto p1 = this->back();
    auto p1_distance = event.distance_from(p1);
    for (auto& p2 : *this) {
      auto p2_distance = event.distance_from(p2);
      if (p1_distance * p2_distance <= 0)
        return true;
      p1_distance = p2_distance;
    }
    return false;
  }

  F max_distance() {
    F distance = 0;
    for (auto& p : *this) {
      distance = compat::max(distance, p.distance_from_origin());
    }
    return distance;
  }

  _ Intersections intersections(Event& event) const {
    auto p1 = this->back();
    auto p1_distance = event.distance_from(p1);
    Intersections intersections;
    for (auto& p2 : *this) {
      auto p2_distance = event.distance_from(p2);
      if (p2_distance == 0) {
        // v2 is crossing point
        intersections.push_back(p2);
      }
      // points lie on opposite sides
      else if (p1_distance * p2_distance < 0) {
        // calculate intersection
        auto direction = p1 - p2;
        auto dot = event.normal.dot(direction);
        auto perpendicular_dot = p1.x * p2.y - p1.y * p2.x;
        intersections.emplace_back((event.distance_from_origin * direction.x -
                                    event.normal.y * perpendicular_dot) /
                                       dot,
                                   (event.distance_from_origin * direction.y +
                                    event.normal.x * perpendicular_dot) /
                                       dot);
      }
      if (intersections.full())
        return intersections;
      p1_distance = p2_distance;
      p1 = p2;
    }
    return intersections;
  }

#if !__CUDACC__

  bool approx_equal(const Polygon& rhs, F epsilon = 1e-5) {
    auto rit = rhs.begin();
    for (auto& p : *this) {
      if (!p.approx_equal(*rit, epsilon))
        return false;
      rit++;
    }
    return true;
  }

  bool approx_equal_circular(const Polygon& rhs, F epsilon = 1e-5) {
    Polygon p(rhs);
    for (int i = 0; i < rhs.size(); i++) {
      std::rotate(p.begin(), p.begin() + 1, p.end());
      if (approx_equal(p, epsilon))
        return true;
    }
    return false;
  }

  bool approx_equal_dihedral(const Polygon& rhs, F epsilon = 1e-5) {
    if (approx_equal_circular(rhs, epsilon))
      return true;
    Polygon r(rhs);

    std::reverse(r.begin(), r.end());

    if (approx_equal_circular(r, epsilon))
      return true;

    return false;
  }

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
