#include "util/test.h"

#include "circle_detector.h"

#include "common/types.h"

using Point = PET2D::Point<F>;
using Event = PET2D::Barrel::Event<F>;
using CircleDetector = PET2D::Barrel::CircleDetector<F>;

TEST("2d/barrel/circle_detector/ctor") {
  CircleDetector circle(0.01);

  CHECK(circle.center.x == 0);
  CHECK(circle.center.y == 0);
}

TEST("2d/barrel/circle_detector/json") {
  CircleDetector circle(1, { 1, 2 });

  json j(circle);
  REQUIRE(j.dump() == "{\"Circle\":{\"center\":[1,2],\"radius\":1}}");
}

TEST("2d/barrel/circle_detector/move") {

  CircleDetector circle(0.01);

  CircleDetector::Vector v(0.5, 0.7);
  circle.center += v;

  CHECK(circle.center.x == F(0.5));
  CHECK(circle.center.y == F(0.7));
  auto phi = M_PI / 6.0;

  CircleDetector rcircle = circle;
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
  CircleDetector circle(1, Point(1, 1));

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
