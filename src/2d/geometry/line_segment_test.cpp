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
  REQUIRE(-distance_from == Approx(std::sqrt(F(5.0))).epsilon(1.0e-7));

  auto t = segment.projection_scaled(p);
  REQUIRE(t == Approx(1 / 3.0).epsilon(1.0e-7));

  auto s = segment.projection(p);
  REQUIRE(s == Approx(std::sqrt(5.0)).epsilon(1.0e-7));

  auto r = segment.projection_relative_middle(p);
  REQUIRE(
      r ==
      Approx(std::sqrt(5.0) - 0.5 * std::sqrt(3 * 3 + 6 * 6)).epsilon(1.0e-7));
}

TEST("2d/geometry/distance_sign") {
  // pointing up
  {
    Point start(0, 0);
    Point end(0, 1);
    LineSegment segment(start, end);

    // point on right
    REQUIRE(-1 == segment.distance_from(Point(1, 0)));
    // point on left
    REQUIRE(1 == segment.distance_from(Point(-1, 0)));
  }
}

TEST("2d/geometry/line_segment/write_read") {
  Point start(3, 0);
  Point end(0, 6);

  LineSegment segment(start, end);
  {
    auto fn = "segment_test.txt"_temp;
    std::ofstream out(fn);
    out << segment;
    out.close();

    std::ifstream in(fn);
    REQUIRE(in);
    LineSegment segment_read(in);
    in.close();

    CHECK(segment.start == segment_read.start);

    REQUIRE(std::remove(fn.c_str()) == 0);
  }
  {
    auto fn = "segment_test.bin"_temp;
    util::obstream out(fn);
    out << segment;
    out.close();

    util::ibstream in(fn);
    REQUIRE(in);
    LineSegment segment_read(in);
    in.close();

    CHECK(segment.start == segment_read.start);

    REQUIRE(std::remove(fn.c_str()) == 0);
  }
}
