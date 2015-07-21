#include "util/test.h"

#include "line_segment.h"

#include "common/types.h"

using Point = PET2D::Point<F>;
using LineSegment = PET2D::LineSegment<F>;

TEST("2d/geometry/line_segment") {
  Point start(3, 0);
  Point end(0, 6);

  LineSegment segment(start, end);

  Point p(4, 3);

  auto distance_from = segment.distance_from(p);
  REQUIRE(distance_from == -Approx(std::sqrt(5.0)).epsilon(1.0e-7));

  auto t = segment.projection_scaled(p);
  REQUIRE(t == Approx(1 / 3.0).epsilon(1.0e-7));

  auto s = segment.projection(p);
  REQUIRE(s == Approx(std::sqrt(5.0)).epsilon(1.0e-7));

  auto r = segment.projection_relative_middle(p);
  REQUIRE(
      r ==
      Approx(std::sqrt(5.0) - 0.5 * std::sqrt(3 * 3 + 6 * 6)).epsilon(1.0e-7));
}
