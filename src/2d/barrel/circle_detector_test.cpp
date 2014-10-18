#include "catch.hpp"

#include "circle_detector.h"

using namespace PET2D;
using namespace PET2D::Barrel;

TEST_CASE("2d/barrel/circle_detector/ctor") {

  CircleDetector<> circle(0.01);

  CHECK(circle.center.x == 0);
  CHECK(circle.center.y == 0);
}

TEST_CASE("2d/barrel/circle_detector/move") {

  CircleDetector<> circle(0.01);

  CircleDetector<>::Point p(0.5, 0.7);
  circle += p;

  CHECK(circle.center.x == 0.5);
  CHECK(circle.center.y == 0.7);
  auto phi = M_PI / 6.0;

  CircleDetector<> rcircle = circle;
  rcircle.rotate(phi);

  auto x = circle.center.x;
  auto y = circle.center.y;
  auto s = std::sin(phi);
  auto c = std::cos(phi);
  auto rx = x * c - y * s;
  auto ry = x * s + y * c;

  CHECK(rx == Approx(rcircle.center.x));
  CHECK(ry == Approx(rcircle.center.y));
}
