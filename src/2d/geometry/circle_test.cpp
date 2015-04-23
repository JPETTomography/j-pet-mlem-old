#include <iostream>
#include <fstream>

#include "util/test.h"

#include "circle.h"

using namespace PET2D;

TEST("2d/geometry/circle/init", "circle initialization") {
  Circle<double> c1(1.);

  CHECK(c1.radius == 1.);
  CHECK(c1.radius2 == 1.);

  Circle<double> c2(std::sqrt(2.));

  CHECK(c2.radius == std::sqrt(2.));  // exact!
  CHECK(c2.radius2 == 2._e13);
}

TEST("2d/geometry/circle/secant") {
  Circle<double> c(1);

  SECTION("angle-0", "0 degrees from (0, 0)") {
    Circle<double>::Event zero(0., 0., 0.);
    auto s = c.secant(zero);

    CHECK(std::min(s[0].x, s[1].x) == -1._e13);
    CHECK(std::max(s[0].x, s[1].x) == 1._e13);

    CHECK(s[0].y == 0.);
    CHECK(s[1].y == 0.);

    auto a = c.secant_angles(zero);
    if (a[0] == Approx(-M_PI))
      a[0] += 2. * M_PI;
    if (a[1] == Approx(-M_PI))
      a[1] += 2. * M_PI;

    CHECK(std::min(a[0], a[1]) == 0._e13);
    CHECK(std::max(a[0], a[1]) == Approx(M_PI));

    Circle<double>::Event one(1., 1., 0.);
    auto sone = c.secant(one);
    CHECK(sone.size() == 1);
    CHECK(sone[0].x == 0.);
    CHECK(sone[0].y == 1.);

    Circle<double>::Event two(2., 2., 0.);
    auto stwo = c.secant(two);
    CHECK(stwo.size() == 0);
  }

  SECTION("angle-90", "90 degrees from (0, 0)") {
    Circle<double>::Event zero90(0., 0., M_PI_2);
    auto s = c.secant(zero90);

    CHECK(s[0].x == 0._e13);
    CHECK(s[1].x == 0._e13);

    CHECK(std::min(s[0].y, s[1].y) == -1.);
    CHECK(std::max(s[0].y, s[1].y) == 1.);

    auto a = c.secant_angles(zero90);
    if (a[0] == Approx(-M_PI))
      a[0] += M_2_PI;
    if (a[1] == Approx(-M_PI))
      a[1] += M_2_PI;

    CHECK(std::min(a[0], a[1]) == Approx(-M_PI_2));
    CHECK(std::max(a[0], a[1]) == Approx(M_PI_2));
  }

  SECTION("angle-45", "45 degrees from (1, 0)") {
    Circle<double>::Event xone45(1., 0., M_PI_4);
    auto s = c.secant(xone45);

    CHECK(std::min(s[0].x, s[1].x) == 0._e13);
    CHECK(std::max(s[0].x, s[1].x) == Approx(xone45.x));

    CHECK(std::min(s[0].y, s[1].y) == -1._e13);
    CHECK(std::max(s[0].y, s[1].y) == Approx(xone45.y));
  }
}

TEST("2d/geometry/circle/secant/math") {
  std::ifstream in("math/secant_test.tab");

  if (!in) {
    WARN("cannot open file `math/secant_test.tab'");
    return;
  }

  double r;
  int n_detectors = 0;
  in >> r;
  in >> n_detectors;

  int line = 0;
  while (true) {

    double x, y, angle;
    double x1, y1;
    double x2, y2;

    in >> x >> y >> angle >> x1 >> y1 >> x2 >> y2;

    if (in.eof())
      break;

    line++;
    Circle<double> c(r);
    Circle<double>::Event event(x, y, angle);

    auto secant = c.secant(event);

    REQUIRE(std::min(secant[0].x, secant[1].x) == Approx(std::min(x1, x2)));
    //    CHECK(std::max(secant[0].x, secant[1].x) == Approx(std::max(x1, x2)));

    //    CHECK(std::min(secant[0].y, secant[1].y) == Approx(std::min(y1, y2)));
    //    CHECK(std::max(secant[0].y, secant[1].y) == Approx(std::max(y1, y2)));

    double angle1, angle2;
    in >> angle1 >> angle2;

    int section1, section2;
    in >> section1 >> section2;

  }
}
