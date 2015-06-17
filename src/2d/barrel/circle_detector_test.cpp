#include "util/test.h"

#include "circle_detector.h"

using namespace PET2D;
using namespace PET2D::Barrel;

TEST("2d/barrel/circle_detector/ctor") {

  CircleDetector<float> circle(0.01);

  CHECK(circle.center.x == 0);
  CHECK(circle.center.y == 0);
}

TEST("2d/barrel/circle_detector/move") {

  CircleDetector<float> circle(0.01);

  CircleDetector<float>::Vector v(0.5, 0.7);
  circle.center += v;

  CHECK(circle.center.x == 0.5f);
  CHECK(circle.center.y == 0.7f);
  auto phi = M_PI / 6.0;

  CircleDetector<float> rcircle = circle;
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

TEST("2d/barrel/circle_detector/intersection") {
  using Event = PET2D::Barrel::Event<float>;
  CircleDetector<float> circle(1, Point<float>(1, 1));

  CHECK(circle.center.x == 1);
  CHECK(circle.center.y == 1);

  // horizontal
  CHECK(true == circle.intersects(Event(1, 1.999, 0)));
  CHECK(true == circle.intersects(Event(9999, 1.999, 0)));
  CHECK(false == circle.intersects(Event(1, 2.001, 0)));
  CHECK(false == circle.intersects(Event(9999, 2.001, 0)));

  // vertical
  CHECK(true == circle.intersects(Event(1.999, 1, M_PI_2)));
  CHECK(true == circle.intersects(Event(1.999, 9999, M_PI_2)));
  CHECK(false == circle.intersects(Event(2.001, 1, M_PI_2)));
  CHECK(false == circle.intersects(Event(2.001, 9999, M_PI_2)));
}
