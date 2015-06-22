#if !defined(_MSC_VER) && !defined(__ICC)
// Disable this test for MSVC & ICC:
// 1. MSVC does not support array initializers
// 2. ICC complains about contructor inheritance

#include <cmath>

#include "util/test.h"

#include "geometry.h"

using Point = Geometry::Point<2, float>;

TEST("geometry/point/2d") {
  REQUIRE(Point(1., 0.).x == Approx(Point(0., 1.).rotate(-M_PI_2)[0]));
  REQUIRE(Point(1., 0.)[1] == Approx(Point(0., 1.).rotate(-M_PI_2)[1]));
}

#endif
