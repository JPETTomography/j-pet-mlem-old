#include "catch.hpp"

#include "circle_detector.h"

TEST_CASE("circle_detector/ctor", "[ctor]") {

  CircleDetector<> circle(0.010,0.010);

  CHECK(circle.center.x == 0);
  CHECK(circle.center.y == 0);
};


TEST_CASE("circle_detector/move", "[ctor]") {

  CircleDetector<> circle(0.010,0.010);

  CircleDetector<>::Point p(0.5,0.7);
  circle+=p;

  CHECK(circle.center.x == 0.5);
  CHECK(circle.center.y == 0.7);
  float phi=M_PI/6.0;
  CircleDetector<> rcircle=circle.rotated(phi);
  float x=circle.center.x;
  float y=circle.center.y;
  float s=sin(phi);
  float c=cos(phi);
  float rx = x*c-y*s;
  float ry = x*s+y*c;

  CHECK(rx==Approx(rcircle.center.x) );
  CHECK(ry==Approx(rcircle.center.y) );

};
