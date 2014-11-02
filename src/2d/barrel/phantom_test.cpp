#include "util/test.h"

#include <cmath>

#include "phantom.h"

using namespace PET2D;
using namespace PET2D::Barrel;

TEST_CASE("2d/barrel/phantom/phantom") {

  Phantom<> phantom;
  phantom.emplace_back(Point<>(1, 1), 2, 2, 0, 0.5);
  phantom.emplace_back(Point<>(0, 1), 1, 0.5, M_PI / 3, 0.75);

  SECTION("intensity") {
    REQUIRE(0.50 == phantom.intensity(Point<>(1, 1)));
    REQUIRE(0.50 == phantom.intensity(Point<>(1.563, -0.8545)));
    REQUIRE(0.00 == phantom.intensity(Point<>(-0.677, -2.5)));

    REQUIRE(0.75 == phantom.intensity(Point<>(-0.328, 0.26)));
    REQUIRE(0.75 == phantom.intensity(Point<>(0.4371, 1.792)));
  }

  SECTION("emit") {
    REQUIRE(false == /**/ phantom.test_emit(Point<>(1, 1), .75));
    REQUIRE(true == /***/ phantom.test_emit(Point<>(1, 1), .45));
    REQUIRE(false == /**/ phantom.test_emit(Point<>(1.563, -0.8545), .51));
    REQUIRE(true == /***/ phantom.test_emit(Point<>(1.563, -0.8545), .10));
    REQUIRE(false == /**/ phantom.test_emit(Point<>(-0.677, -2.5), .1));
    REQUIRE(false == /**/ phantom.test_emit(Point<>(-0.677, -2.5), .25));
    REQUIRE(false == /**/ phantom.test_emit(Point<>(-0.677, -2.5), 0.001));

    REQUIRE(false == /**/ phantom.test_emit(Point<>(-0.328, 0.26), .76));
    REQUIRE(true == /***/ phantom.test_emit(Point<>(-0.328, 0.26), .74));
    REQUIRE(false == /**/ phantom.test_emit(Point<>(0.4371, 1.792), .77));
    REQUIRE(true == /***/ phantom.test_emit(Point<>(0.4371, 1.792), .73));
  }
}
