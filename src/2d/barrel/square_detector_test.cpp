#include "util/test.h"

#include "square_detector.h"

using namespace PET2D;
using namespace PET2D::Barrel;

TEST("2d/barrel/square_detector/intersection") {

  SquareDetector<> d(2., 1., 0.);

  CHECK(d[0].x == 1.);
  CHECK(d[0].y == .5);
  CHECK(d[1].x == 1.);
  CHECK(d[1].y == -.5);
  CHECK(d[2].x == -1.);
  CHECK(d[2].y == -.5);
  CHECK(d[3].x == -1.);
  CHECK(d[3].y == .5);

  SquareDetector<>::Event e1(1., 0., M_PI_4);
  SquareDetector<>::Event e2(1., -3., -M_PI_4);

  CHECK(true == d.intersects(e1));
  CHECK(false == d.intersects(e2));

  SECTION("detector/intersection/points", "intersection points") {
    auto i1 = d.intersections(e1);

    REQUIRE(i1.size() == 2);
    CHECK(std::min(i1[0].x, i1[1].x) == 0.5_e13);
    CHECK(std::max(i1[0].x, i1[1].x) == 1.0_e13);

    CHECK(std::min(i1[0].y, i1[1].y) == -0.5_e13);
    CHECK(std::max(i1[0].y, i1[1].y) == 0.0_e13);
  }

  SECTION("detector/rotated", "rotated") {
    auto dr = d;
    dr.rotate(M_PI_4);
    auto s = std::sin(M_PI_4);
    auto c = std::sin(M_PI_4);

    CHECK(dr[0].x == Approx(d[0].x * c - d[0].y * s));
    CHECK(dr[0].y == Approx(d[0].x * s + d[0].y * c));
  }
#if DONT_TEST
  SECTION("detector/translated+rotated", "translated and rotated") {
    Detector<>::Point p(2., 3.);
    auto dtr = (d + p);
    dtr.rotate(M_PI_4);
    auto s = std::sin(M_PI_4);
    auto c = std::sin(M_PI_4);

    CHECK(dtr[0].x == Approx((d[0].x + p.x) * c - (d[0].y + p.y) * s));
    CHECK(dtr[0].y == Approx((d[0].x + p.x) * s + (d[0].y + p.y) * c));
  }
#endif
}
