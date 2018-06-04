#include "util/test.h"

#include <cmath>

#include "event.h"

using Point = PET2D::Point<double>;
using Event = PET2D::Barrel::Event<double>;

TEST("2d/barrel/event/set") {

  Event event(1.0, 0.5, 2.0);

  REQUIRE(1.0 == event.origin.x);
  REQUIRE(0.5 == event.origin.y);
  REQUIRE(event.direction.x == std::cos(2.0));
}

TEST("2d/barrel/event/distance") {
  // pointing up
  {
    PET2D::Barrel::Event<double> event(2, 1, M_PI_2);  // 90 deg
    REQUIRE(2 == event.origin.x);
    REQUIRE(1 == event.origin.y);
    // point on left
    REQUIRE(-1 == event.distance_from(Point(1, 1)));
    // point on right
    REQUIRE(1 == event.distance_from(Point(3, 1)));
  }
  // pointing right
  {
    PET2D::Barrel::Event<double> event(2, 1, 0);
    REQUIRE(2 == event.origin.x);
    REQUIRE(1 == event.origin.y);
    // point above (left)
    REQUIRE(-1 == event.distance_from(Point(1, 2)));
    // point below (right)
    REQUIRE(1 == event.distance_from(Point(1, 0)));
  }
}
