#if !defined(_MSC_VER) && !defined(__ICC)
// Disable this test for MSVC & ICC:
// 1. MSVC does not support array initializers
// 2. ICC complains about contructor inheritance

#include <cmath>

#include "util/test.h"

#include "geometry.h"

using namespace Geometry;

TEST_CASE("geometry/point/2d") {
  REQUIRE(Point<2>(1., 0.) == Point<2>(1., 0.));
  REQUIRE(Point<2>(1., 0.)[0] == Approx(Point<2>(0., 1.).rotate(-M_PI_2)[0]));
  REQUIRE(Point<2>(1., 0.)[1] == Approx(Point<2>(0., 1.).rotate(-M_PI_2)[1]));
}

#endif
