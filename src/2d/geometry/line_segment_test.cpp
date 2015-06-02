#include "util/test.h"

#include "line_segment.h"

TEST("PET2D/geometry/line_segment") {
  using Point = PET2D::Point<float>;

  Point start(3, 0);
  Point end(0, 6);

  PET2D::LineSegment<float> segment(start, end);

  Point p(4, 3);

  auto distance_from = segment.distance_from(p);
  REQUIRE(distance_from == -Approx(std::sqrt(5.0)).epsilon(1.0e-7));

  auto t = segment.projection_scaled(p);
  REQUIRE(t == Approx(1 / 3.0).epsilon(1.0e-7));
}
