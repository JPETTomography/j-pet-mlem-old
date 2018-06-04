#pragma once

#if !__CUDACC__
#include "util/bstream.h"
#endif

#include "2d/geometry/vector.h"
#include "2d/geometry/point.h"

namespace PET2D {

/// Line segment
////
/// Gives some extra functions allowing calculation Point distance to the
/// segment. The distance is calculated using normalized normal vector
/// \f$ (A, B) \f$ and distance to origin \f$ C \f$ using formula:
/// \f[
///   A x + B y + C
/// \f]
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

#if !__CUDACC__
  LineSegment(std::istream& in) : LineSegment(Point(in), in) {}
  LineSegment(util::ibstream& in) : LineSegment(Point(in), in) {}

 private:
  // This private proxy constructors guarantee proper read order:
  LineSegment(const Point&& start, std::istream& in)
      : LineSegment(start, Point(in)) {}
  LineSegment(const Point&& start, util::ibstream& in)
      : LineSegment(start, Point(in)) {}

 public:
  friend std::ostream& operator<<(std::ostream& out,
                                  const LineSegment& segment) {
    out << segment.start << ' ' << segment.end;
    return out;
  }
  friend util::obstream& operator<<(util::obstream& out,
                                    const LineSegment& segment) {
    out << segment.start << segment.end;
    return out;
  }
#endif

  /// Returns distance from given Point to the line containing LineSegment.
  ////
  /// \return 0 means points lies on the line, 1 left, -1 right
  /// \todo FIXME: This is different from Barrel::Event::distance_from returning
  /// opposite.
  _ F distance_from(const Point& p) const {
    return p.as_vector().dot(normal) - distance_from_origin;
  }

  /// Returns Point projected to the line containing LineSegment.
  ////
  /// Projection space origin is at the start point of the segment.
  _ F projection(const Point& p) const { return (p - start).dot(direction); }

  /// Returns Point projected to the line containing LineSegment, length scaled.
  ////
  /// Projection space origin is at the start point of the segment.
  _ F projection_scaled(const Point& p) const { return projection(p) / length; }

  /// Returns Point projected to the line containing LineSegment.
  ////
  /// Projection space origin is at the mid point of the segment.
  _ F projection_relative_middle(const Point& p) const {
    return (p - mid_point).dot(direction);
  }

  // FIXME: these should be const, but making them const breaks too much code
  Point start;             ///< line segment 1st endpoint
  Point end;               ///< line segment 2nd endpoint
  Point mid_point;         ///< point lying exactly in the middle of the segment
  Vector direction;        ///< direction vector from start towards end
  Vector normal;           ///< vector normal to the direction
  F length;                ///< length of the segment
  F distance_from_origin;  ///< distance of the segment from the origin
};

}  // PET2D
