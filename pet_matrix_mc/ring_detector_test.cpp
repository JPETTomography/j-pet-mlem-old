#include <catch.hpp>

#include <cmath>

#include "ring_detector.h"

TEST_CASE("ring_detector/center", "ring detector center") {
  event<> event(0.0, 0.0, 2.0);

  auto time = tof(event, 1.0);

  REQUIRE(  1.0 == time.first  );
  REQUIRE( -1.0 == time.second );
}

TEST_CASE("ring_detector/on_x_axis", "ring detector x axis") {
  event<> event1(0.5, 0.0, 0.0);

  auto time = tof(event1, 1.0);

  REQUIRE(  0.5 == time.first  );
  REQUIRE( -1.5 == time.second );

  event<> event2(0.5, 0.0, M_PI);

  time = tof(event2, 1.0);

  REQUIRE(  1.5 == time.first  );
  REQUIRE( -0.5 == time.second );
}

TEST_CASE("ring_detector/on_y_axis", "ring detector y axis") {
  event<> event1(0.0, 0.5, M_PI/2.0);

  auto time = tof(event1, 1.0);

  REQUIRE(  0.5 == time.first  );
  REQUIRE( -1.5 == time.second );

  event<> event2(0.0, 0.5, -M_PI/2.0);

  time = tof(event2, 1.0);

  REQUIRE(  1.5 == time.first  );
  REQUIRE( -0.5 == time.second );
}

TEST_CASE("ring_detector/on_sym", "ring detector symmetry") {
  auto x = 0.5;
  auto y = 0.5;

  event<> event1(x, y, M_PI_4);

  auto time = tof(event1, 1.0);

  auto dist = std::sqrt(x*x+y*y);

  REQUIRE( Approx( 1.0 - dist) == time.first  );
  REQUIRE( Approx(-1.0 - dist) == time.second );

  event<> event2(x, y, 5.0*M_PI_4);

  time = tof(event2, 1.0);

  REQUIRE( Approx( 1.0 + dist) == time.first  );
  REQUIRE( Approx(-1.0 + dist) == time.second );
}

TEST_CASE("ring_detector/radius", "ring detector radious") {
  auto x =  0.8;
  auto y = -0.6;

  event<> event(x, y,2.1);

  auto time = tof(event, 1.0);

  auto s = std::sin(event.phi);
  auto c = std::cos(event.phi);

  auto x_hit = x+time.first*c;
  auto y_hit = y+time.first*s;

  REQUIRE( 1.0 == Approx(x_hit*x_hit+y_hit*y_hit) );

  x_hit = x+time.second*c;
  y_hit = y+time.second*s;

  REQUIRE( 1.0 == Approx(x_hit*x_hit+y_hit*y_hit) );
}

TEST_CASE("ring_detector/lor_center", "LOR center") {
  event<> event1(0.0, 0.0, 0.001);

  auto time = tof(event1, 1.0);
  auto lors = lor(time, event1, 1.0, 128);

  REQUIRE(  0 == lors.first  );
  REQUIRE( 64 == lors.second );

  event<> event2(0.0, 0.0, M_PI/2.0+0.001);

  time = tof(event2, 1.0);
  lors = lor(time, event2, 1.0, 128);

  REQUIRE( 32 == lors.first  );
  REQUIRE( 96 == lors.second );

  event<> event3(0.0, 0.0, 3.0*M_PI/2.0+0.001);

  time = tof(event3, 1.0);
  lors = lor(time, event3, 1.0, 128);

  REQUIRE( 32 == lors.first  );
  REQUIRE( 96 == lors.second );
}
