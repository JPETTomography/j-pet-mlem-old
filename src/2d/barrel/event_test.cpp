#include "util/test.h"

#include <cmath>

#include "event.h"

using namespace PET2D;
using namespace PET2D::Barrel;

TEST("2d/barrel/event/set") {

  Event<double> event(1.0, 0.5, 2.0);

  REQUIRE(1.0 == event.x);
  REQUIRE(0.5 == event.y);
  REQUIRE(event.direction.x == std::cos(2.0));
}


