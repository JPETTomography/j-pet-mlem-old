#include "catch.hpp"

#include <cmath>

#include "phantom.h"

TEST_CASE("2d_xy/phantom/elliptical_region") {

  EllipticalRegion<> disk(Point<>(1, 1), 2, 2, 0, 0.5);
  EllipticalRegion<> region(Point<>(0, 1), 1, 0.5, M_PI / 3, 0.75);

  SECTION("getter") {
    REQUIRE(disk.activity() == 0.5);
    REQUIRE(region.activity() == 0.75);
  }

  SECTION("operator()") {
    REQUIRE(true == /***/ disk(Point<>(1, 1)));
    REQUIRE(true == /***/ disk(Point<>(1.563, -0.8545)));
    REQUIRE(false == /**/ disk(Point<>(-0.677, -2.5)));

    REQUIRE(true == /***/ region(Point<>(-0.328, 0.26)));
    REQUIRE(true == /***/ region(Point<>(0.4371, 1.792)));
    REQUIRE(false == /**/ region(Point<>(1, 1)));
  }
}

TEST_CASE("2d_xy/phantom/phantom") {

  Phantom<> phantom;
  phantom.push_back(EllipticalRegion<>(Point<>(1, 1), 2, 2, 0, 0.5));
  phantom.push_back(EllipticalRegion<>(Point<>(0, 1), 1, 0.5, M_PI / 3, 0.75));

  SECTION("activity") {
    REQUIRE(0.50 == phantom.activity(Point<>(1, 1)));
    REQUIRE(0.50 == phantom.activity(Point<>(1.563, -0.8545)));
    REQUIRE(0.00 == phantom.activity(Point<>(-0.677, -2.5)));

    REQUIRE(0.75 == phantom.activity(Point<>(-0.328, 0.26)));
    REQUIRE(0.75 == phantom.activity(Point<>(0.4371, 1.792)));
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
