#include "catch.hpp"

#include <cmath>

#include "phantom.h"

TEST_CASE("elliptical_region", "elliptical region") {

  EllipticalRegion disk(1.0, 1.0, 2.0, 2.0, 0.0, 0.5);
  EllipticalRegion region(0, 1, 1, 0.5, M_PI / 3.0, 0.75);

  SECTION("getter", "") {
    REQUIRE(disk.activity() == 0.5);
    REQUIRE(region.activity() == 0.75);
  }

  SECTION("in", "") {
    REQUIRE(true == disk.in(1, 1));
    REQUIRE(true == disk.in(1.563, -0.8545));
    REQUIRE(false == disk.in(-0.677, -2.5));

    REQUIRE(true == region.in(-0.328, 0.26));
    REQUIRE(true == region.in(0.4371, 1.792));
    REQUIRE(false == region.in(1, 1));
  }
}

TEST_CASE("phantom", "phantom") {

  // Phantom phantom(-3, -3, 3, 3);
  Phantom phantom;
  phantom.addRegion(1.0, 1.0, 2.0, 2.0, 0.0, 0.5);
  phantom.addRegion(0, 1, 1, 0.5, M_PI / 3.0, 0.75);

  SECTION("activity", "") {
    REQUIRE(0.5 == phantom.activity(1, 1));
    REQUIRE(0.5 == phantom.activity(1.563, -0.8545));
    REQUIRE(0.0 == phantom.activity(-0.677, -2.5));

    REQUIRE(0.75 == phantom.activity(-0.328, 0.26));
    REQUIRE(0.75 == phantom.activity(0.4371, 1.792));
  }

  SECTION("emit", "") {
    REQUIRE(false == phantom.emit(1, 1, .75));
    REQUIRE(true == phantom.emit(1, 1, .45));
    REQUIRE(false == phantom.emit(1.563, -0.8545, .51));
    REQUIRE(true == phantom.emit(1.563, -0.8545, .10));
    REQUIRE(false == phantom.emit(-0.677, -2.5, .1));
    REQUIRE(false == phantom.emit(-0.677, -2.5, .25));
    REQUIRE(false == phantom.emit(-0.677, -2.5, 0.001));

    REQUIRE(false == phantom.emit(-0.328, 0.26, .76));
    REQUIRE(true == phantom.emit(-0.328, 0.26, .74));
    REQUIRE(false == phantom.emit(0.4371, 1.792, .77));
    REQUIRE(true == phantom.emit(0.4371, 1.792, .73));
  }
}
