#include "util/test.h"

#include "rectangle.h"

#include "common/types.h"

using Point = PET2D::Point<F>;
using Rectangle = PET2D::Rectangle<F>;
using RectanglePointGenerator = PET2D::RectanglePointGenerator<F>;

TEST("2d/geometry/rectangle") {
  Rectangle r(1, 2, 3, 4);

  REQUIRE(r.area == 48.0_e7);
  REQUIRE(r.contains(Point(1, 2)));
  REQUIRE(r.contains(Point(1, 5.9)));
  REQUIRE(r.contains(Point(-1.99, 2)));

  REQUIRE_FALSE(r.contains(Point(1, 6.1)));
  REQUIRE_FALSE(r.contains(Point(-2.1, 2)));
}

TEST("2d/geometry/rectangle/point_generator") {
  Rectangle r(1, 2, 3, 4);
  RectanglePointGenerator gen(r);
  std::mt19937_64 rng;

  for (int i = 0; i < 100; i++) {
    auto p = gen(rng);
    REQUIRE(r.contains(p));
  }
}
