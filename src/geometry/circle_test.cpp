#include <iostream>
#include <fstream>

#include "catch.hpp"

#include "circle.h"

TEST_CASE("circle/init", "circle initialization") {
  Circle<> c1(1.);

  CHECK(c1.radius() == 1.);
  CHECK(c1.radius2() == 1.);

  Circle<> c2(std::sqrt(2.));

  CHECK(c2.radius() == std::sqrt(2.));  // exact!
  CHECK(c2.radius2() == Approx(2.));
}

TEST_CASE("circle/secant") {
  Circle<> c(1);

  SECTION("angle-0", "0 degrees from (0, 0)") {
    Circle<>::Event zero(0., 0., 0.);
    auto s = c.secant(zero);

    CHECK(std::min(s[0].x, s[1].x) == Approx(-1.));
    CHECK(std::max(s[0].x, s[1].x) == Approx(1.));

    CHECK(s[0].y == 0.);
    CHECK(s[1].y == 0.);

    auto a = c.secant_angles(zero);
    if (a[0] == Approx(-M_PI))
      a[0] += 2. * M_PI;
    if (a[1] == Approx(-M_PI))
      a[1] += 2. * M_PI;

    CHECK(std::min(a[0], a[1]) == Approx(0.));
    CHECK(std::max(a[0], a[1]) == Approx(M_PI));
  }
  SECTION("angle-90", "90 degrees from (0, 0)") {
    Circle<>::Event zero90(0., 0., M_PI_2);
    auto s = c.secant(zero90);

    CHECK(s[0].x == Approx(0.));
    CHECK(s[1].x == Approx(0.));

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
    Circle<>::Event xone45(1., 0., M_PI_4);
    auto s = c.secant(xone45);

    CHECK(std::min(s[0].x, s[1].x) == Approx(0.0).epsilon(1.0e-13));
    CHECK(std::max(s[0].x, s[1].x) == Approx(xone45.x));

    CHECK(std::min(s[0].y, s[1].y) == Approx(-1.0));
    CHECK(std::max(s[0].y, s[1].y) == Approx(xone45.y));
  }
}

TEST_CASE("circle/secant/math", "[math]") {
  std::ifstream in("secant.test");

  if (!in) {
    WARN("cannot open file `secant.test'");
    return;
  }

  double r;
  int n_detectors;
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
    Circle<>::Event event(x, y, angle);

    auto secant = c.secant(event);

    CHECK(std::min(secant[0].x, secant[1].x) == Approx(std::min(x1, x2)));
    CHECK(std::max(secant[0].x, secant[1].x) == Approx(std::max(x1, x2)));

    CHECK(std::min(secant[0].y, secant[1].y) == Approx(std::min(y1, y2)));
    CHECK(std::max(secant[0].y, secant[1].y) == Approx(std::max(y1, y2)));

    double angle1, angle2;
    in >> angle1 >> angle2;
    auto angles = c.secant_angles(event);

    CHECK(std::min(angles[0], angles[1]) == Approx(std::min(angle1, angle2)));
    CHECK(std::max(angles[0], angles[1]) == Approx(std::max(angle1, angle2)));

    int section1, section2;
    in >> section1 >> section2;
    auto sections = c.secant_sections(event, n_detectors);

    CHECK(std::min(sections[0], sections[1]) == std::min(section1, section2));
    CHECK(std::max(sections[0], sections[1]) == std::max(section1, section2));
  }
}