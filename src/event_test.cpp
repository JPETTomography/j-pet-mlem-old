#include <catch.hpp>

#include <cmath>

#include "event.h"

TEST_CASE("event/set", "event set") {
  Event<double> event(1.0, 0.5, 2.0);

  REQUIRE(1.0 == event.x);
  REQUIRE(0.5 == event.y);
  REQUIRE(2.0 == event.phi);
}
