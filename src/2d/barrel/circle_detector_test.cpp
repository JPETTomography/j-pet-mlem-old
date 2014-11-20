#include "util/test.h"

#include "circle_detector.h"

using namespace PET2D;
using namespace PET2D::Barrel;

TEST("2d/barrel/circle_detector/ctor") {

  CircleDetector<> circle(0.01);

  CHECK(circle.x == 0);
  CHECK(circle.y == 0);
}

TEST("2d/barrel/circle_detector/move") {

  CircleDetector<> circle(0.01);

  CircleDetector<>::Point p(0.5, 0.7);
  circle += p;

  CHECK(circle.x == 0.5);
  CHECK(circle.y == 0.7);
  auto phi = M_PI / 6.0;

  CircleDetector<> rcircle = circle;
  rcircle.rotate(phi);

  auto x = circle.x;
  auto y = circle.y;
  auto s = std::sin(phi);
  auto c = std::cos(phi);
  auto rx = x * c - y * s;
  auto ry = x * s + y * c;

  CHECK(rx == Approx(rcircle.x));
  CHECK(ry == Approx(rcircle.y));
}
