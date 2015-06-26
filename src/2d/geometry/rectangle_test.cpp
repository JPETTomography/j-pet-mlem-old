#include "util/test.h"

#include "rectangle.h"

TEST("2d/geometry/rectangle") {
  using Point = PET2D::Point<float>;

  PET2D::Rectangle<float> r(1, 2, 3, 4);

  REQUIRE(r.area == 48.0_e7);
  REQUIRE(r.contains(Point(1, 2)));
  REQUIRE(r.contains(Point(1, 5.9)));
  REQUIRE(r.contains(Point(-1.99, 2)));

  REQUIRE(!r.contains(Point(1, 6.1)));
  REQUIRE(!r.contains(Point(-2.1, 2)));

}
