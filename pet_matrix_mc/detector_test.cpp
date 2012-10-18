#include <catch.hpp>

#include <cmath>

#include "detector.h"

TEST_CASE("detector/center", "detector center") {
  event<double> event(0.0, 0.0, 2.0);

  std::pair<double, double> time = tof(event, 1.0);

  REQUIRE(  1.0 == time.first  );
  REQUIRE( -1.0 == time.second );
}

TEST_CASE("detector/on_x_axis", "detector x axis") {
  event<double> event1(0.5, 0.0, 0.0);

  std::pair<double, double> time = tof(event1, 1.0);

  REQUIRE(  0.5 == time.first  );
  REQUIRE( -1.5 == time.second );

  event<double> event2(0.5, 0.0, M_PI);

  time = tof(event2, 1.0);

  REQUIRE(  1.5 == time.first  );
  REQUIRE( -0.5 == time.second );
}

TEST_CASE("detector/on_y_axis", "detector y axis") {
  event<double> event1(0.0, 0.5, M_PI/2.0);

  std::pair<double, double> time = tof(event1, 1.0);

  REQUIRE(  0.5 == time.first  );
  REQUIRE( -1.5 == time.second );

  event<double> event2(0.0, 0.5, -M_PI/2.0);

  time = tof(event2, 1.0);

  REQUIRE(  1.5 == time.first  );
  REQUIRE( -0.5 == time.second );
}

TEST_CASE("detector/on_sym", "detector symmetry") {
  double x = 0.5;
  double y = 0.5;

  event<double> event1(x, y, M_PI_4);

  std::pair<double, double> time = tof(event1, 1.0);

  double dist = sqrt(x*x+y*y);

  REQUIRE( Approx( 1.0 - dist) == time.first  );
  REQUIRE( Approx(-1.0 - dist) == time.second );

  event<double> event2(x, y, 5.0*M_PI_4);

  time = tof(event2, 1.0);

  REQUIRE( Approx( 1.0 + dist) == time.first  );
  REQUIRE( Approx(-1.0 + dist) == time.second );
}

TEST_CASE("detector/radius", "detector radious") {
  double x= 0.8;
  double y = -0.6;

  event<double> event(x, y,2.1);

  std::pair<double, double> time = tof(event, 1.0);

  double s, c;
  sincos(event.phi(), s, c);

  double x_hit = x+time.first*c;
  double y_hit = y+time.first*s;

  REQUIRE( 1.0 == Approx(x_hit*x_hit+y_hit*y_hit) );

  x_hit = x+time.second*c;
  y_hit = y+time.second*s;

  REQUIRE( 1.0 == Approx(x_hit*x_hit+y_hit*y_hit) );
}

TEST_CASE("detector/lor_center", "LOR center") {
  event<double> event1(0.0, 0.0, 0.001);

  std::pair<double, double> time = tof(event1, 1.0);
  std::pair<short, short> lors = lor(time, event1, 1.0, 128);

  REQUIRE(  0 == lors.first  );
  REQUIRE( 64 == lors.second );

  event<double> event2(0.0, 0.0, M_PI/2.0+0.001);

  time = tof(event2, 1.0);
  lors = lor(time, event2, 1.0, 128);

  REQUIRE( 32 == lors.first  );
  REQUIRE( 96 == lors.second );

  event<double> event3(0.0, 0.0, 3.0*M_PI/2.0+0.001);

  time = tof(event3, 1.0);
  lors = lor(time, event3, 1.0, 128);

  REQUIRE( 32 == lors.first  );
  REQUIRE( 96 == lors.second );
}
