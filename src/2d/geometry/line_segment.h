#ifndef LINE_SEGMENT
#define LINE_SEGMENT

#include "2d/geometry/vector.h"
#include "2d/geometry/point.h"

namespace PET2D {
template <typename FType> class LineSegment {
  using F = FType;
  using Point = PET2D::Point<F>;
  using Vector = PET2D::Vector<F>;

 public:
  LineSegment(const Point& start, const Point& end)
      : start(start),
        end(end),
        mid_point(Point((start.x + end.x) / 2, (start.y + end.y) / 2)),
        direction((end - start).normalized()),
        normal(direction.perpendicular()),
        distance(dot(start.as_vector(), normal)) {}

  F distance_from(const Point& p);
  F projection_on(const Point& p);

  const Point start;
  const Point end;
  const Point mid_point;
  const Vector direction;
  const Vector normal;
  const F distance;
};
}

#endif  // LINE_SEGMENT
