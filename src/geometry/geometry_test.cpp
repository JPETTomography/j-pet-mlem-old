#include <cmath>

#include "catch.hpp"

#include "geometry.h"

using namespace Geometry;

TEST_CASE("geometry/point/2d") {
  REQUIRE(Point<2>(1., 0.) == Point<2>(1., 0.));
  REQUIRE(Point<2>(1., 0.)[0] == Approx(Point<2>(0., 1.).rotate(-M_PI_2)[0]));
  REQUIRE(Point<2>(1., 0.)[1] == Approx(Point<2>(0., 1.).rotate(-M_PI_2)[1]));
}
