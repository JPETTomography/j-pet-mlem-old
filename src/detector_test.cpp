#include "catch.hpp"

#include "detector.h"

TEST_CASE("detector/intersection", "polygon intersection") {
  Detector<> d(2., 1.);

  CHECK(d[0].x == 1.);
  CHECK(d[0].y == .5);
  CHECK(d[1].x == 1.);
  CHECK(d[1].y == -.5);
  CHECK(d[2].x == -1.);
  CHECK(d[2].y == -.5);
  CHECK(d[3].x == -1.);
  CHECK(d[3].y == .5);

  Detector<>::Event e1(1., 0., M_PI_4);
  Detector<>::Event e2(1., -3., -M_PI_4);

  CHECK(true == d.intersects(e1));
  CHECK(false == d.intersects(e2));

  SECTION("detector/intersection/points", "intersection points") {
    auto i1 = d.intersections(e1);

    REQUIRE(i1.size() == 2);
    CHECK(std::min(i1[0].x, i1[1].x) == Approx(0.5));
    CHECK(std::max(i1[0].x, i1[1].x) == Approx(1.));

    CHECK(std::min(i1[0].y, i1[1].y) == Approx(-0.5));
    CHECK(std::max(i1[0].y, i1[1].y) == Approx(0.));
  }

  SECTION("detector/rotated", "rotated") {
    auto dr = d.rotated(M_PI_4);
    auto s = std::sin(M_PI_4);
    auto c = std::sin(M_PI_4);

    CHECK(dr[0].x == Approx(d[0].x * c - d[0].y * s));
    CHECK(dr[0].y == Approx(d[0].x * s + d[0].y * c));
  }
#if 0
  SECTION("detector/translated+rotated", "translated and rotated") {
    Detector<>::Point p(2., 3.);
    auto dtr = (d + p).rotated(M_PI_4);
    auto s = std::sin(M_PI_4);
    auto c = std::sin(M_PI_4);

    CHECK(dtr[0].x == Approx((d[0].x + p.x) * c - (d[0].y + p.y) * s));
    CHECK(dtr[0].y == Approx((d[0].x + p.x) * s + (d[0].y + p.y) * c));
  }
#endif
}
