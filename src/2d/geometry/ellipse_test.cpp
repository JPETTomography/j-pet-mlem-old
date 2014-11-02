#include "util/test.h"

#include <cmath>

#include "ellipse.h"

using namespace PET2D;

TEST("2d/geometry/ellipse/elliptical_region") {

  EllipticalSource<> s1(1, 1, 2, 2, 0, 0.5);
  EllipticalSource<> s2(0, 1, 1, 0.5, M_PI / 3, 0.75);

  SECTION("getter") {
    REQUIRE(s1.intensity == 0.5);
    REQUIRE(s2.intensity == 0.75);
  }

  SECTION("contains") {
    REQUIRE(true == /***/ s1.contains({ 1, 1 }));
    REQUIRE(true == /***/ s1.contains({ 1.563, -0.8545 }));
    REQUIRE(false == /**/ s1.contains({ -0.677, -2.5 }));

    REQUIRE(true == /***/ s2.contains({ -0.328, 0.26 }));
    REQUIRE(true == /***/ s2.contains({ 0.4371, 1.792 }));
    REQUIRE(false == /**/ s2.contains({ 1, 1 }));
  }
}
