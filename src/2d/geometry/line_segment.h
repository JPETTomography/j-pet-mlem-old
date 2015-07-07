#pragma once

#include "2d/geometry/vector.h"
#include "2d/geometry/point.h"

namespace PET2D {

template <typename FType> struct LineSegment {
  using F = FType;
  using Point = PET2D::Point<F>;
  using Vector = PET2D::Vector<F>;

  LineSegment() = default;

  LineSegment(const Point& start, const Point& end)
      : start(start),
        end(end),
        mid_point(Point((start.x + end.x) / 2, (start.y + end.y) / 2)),
        direction((end - start).normalized()),
        normal(direction.perpendicular()),
        length((end - start).length()),
        distance_from_origin(start.as_vector().dot(normal)) {}

  LineSegment(std::istream& in) : LineSegment(Point(in), Point(in)) {}

  F distance_from(const Point& p) {
    return p.as_vector().dot(normal) - distance_from_origin;
  }

  F projection(const Point& p) { return (p - start).dot(direction); }
  F projection_scaled(const Point& p) { return projection(p) / length; }

  F projection_relative_middle(const Point& p) {
    return projection(p) - 0.5 * length;
  }

  Point start;
  Point end;
  Point mid_point;
  Vector direction;
  Vector normal;
  F length;
  F distance_from_origin;
};
}
