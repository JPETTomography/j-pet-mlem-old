#include <iostream>
#include <fstream>

#include "catch.hpp"

#include "polygon.h"

TEST_CASE("polygon/intersection", "polygon intersection") {
  Polygon<> p;
  p.push_back(Point<>(1., 1.));
  p.push_back(Point<>(2., 1.));
  p.push_back(Point<>(2., 2.));
  p.push_back(Point<>(1., 2.));

  Polygon<>::Event e1(0., .5, M_PI_4);
  Polygon<>::Event e2(0., 1.5, M_PI_4);
  Polygon<>::Event e3(0., 1.5, 0.);

  CHECK(true == p.intersects(e1));
  CHECK(false == p.intersects(e2));
  CHECK(true == p.intersects(e3));

  SECTION("polygon/intersection/points", "intersection points") {
    auto i1 = p.intersections(e1);

    REQUIRE(i1.size() == 2);
    CHECK(std::min(i1[0].x, i1[1].x) == Approx(1.));
    CHECK(std::max(i1[0].x, i1[1].x) == Approx(1.5));

    CHECK(std::min(i1[0].y, i1[1].y) == Approx(1.5));
    CHECK(std::max(i1[0].y, i1[1].y) == Approx(2.));

    auto i3 = p.intersections(e3);

    REQUIRE(i3.size() == 2);
    CHECK(std::min(i3[0].x, i3[1].x) == Approx(1.));
    CHECK(std::max(i3[0].x, i3[1].x) == Approx(2.));

    CHECK(std::min(i3[0].y, i3[1].y) == Approx(1.5));
    CHECK(std::max(i3[0].y, i3[1].y) == Approx(1.5));
  }
}

TEST_CASE("polygon/intersection/math", "rectangle inters from mathematica") {
  std::ifstream in("polygon.test");

  if (!in) {
    WARN("cannot open file `polygon.test'");
    return;
  }

  int n_events;
  in >> n_events;

  Polygon<> poly;

  for (int i = 0; i < 4; ++i) {
    double x, y;
    in >> x >> y;
    Polygon<>::Point p(x, y);
    poly.push_back(p);
  }

  for (int i = 0; i < n_events; i++) {
    double x, y, phi;
    in >> x >> y >> phi;

    double a, b, c;
    in >> a >> b >> c;

    size_t n_iters;
    in >> n_iters;

    Polygon<>::Event event(x, y, phi);
    bool intersects = n_iters > 0;

    CHECKED_IF(poly.intersects(event) == intersects) {
      auto inters = poly.intersections(event);

      CHECKED_IF(inters.size() == n_iters) {

        if (n_iters > 0) {
          double ix, iy;
          Polygon<>::Intersections m_inters;

          for (size_t j = 0; j < n_iters; ++j) {
            in >> ix >> iy;
            Polygon<>::Point p(ix, iy);
            m_inters.push_back(p);
          }

          if (n_iters == 1) {
            CHECK(inters[0].x == Approx(m_inters[0].x));
            CHECK(inters[0].y == Approx(m_inters[0].y));
          } else {
            CHECK(std::min(inters[0].x, inters[1].x) ==
                  Approx(std::min(m_inters[0].x, m_inters[1].x)));
            CHECK(std::min(inters[0].y, inters[1].y) ==
                  Approx(std::min(m_inters[0].y, m_inters[1].y)));

            CHECK(std::max(inters[0].x, inters[1].x) ==
                  Approx(std::max(m_inters[0].x, m_inters[1].x)));
            CHECK(std::max(inters[0].y, inters[1].y) ==
                  Approx(std::max(m_inters[0].y, m_inters[1].y)));
          }
        }
      }
    }
  }
}
