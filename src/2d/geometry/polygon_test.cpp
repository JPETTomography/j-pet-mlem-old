#include <iostream>
#include <fstream>

#include "util/test.h"

#include "polygon.h"

using namespace PET2D;

TEST_CASE("geometry/2d/barrel/polygon/intersection") {
  Polygon<4> p;
  p.emplace_back(1., 1.);
  p.emplace_back(2., 1.);
  p.emplace_back(2., 2.);
  p.emplace_back(1., 2.);

  Polygon<4>::Event e1(0., .5, M_PI_4);
  Polygon<4>::Event e2(0., 1.5, M_PI_4);
  Polygon<4>::Event e3(0., 1.5, 0.);

  CHECK(true == p.intersects(e1));
  CHECK(false == p.intersects(e2));
  CHECK(true == p.intersects(e3));

  SECTION("polygon/intersection/points", "intersection points") {
    auto i1 = p.intersections(e1);

    REQUIRE(i1.size() == 2);
    CHECK(std::min(i1[0].x, i1[1].x) == 1.0_e13);
    CHECK(std::max(i1[0].x, i1[1].x) == 1.5_e13);

    CHECK(std::min(i1[0].y, i1[1].y) == 1.5_e13);
    CHECK(std::max(i1[0].y, i1[1].y) == 2.0_e13);

    auto i3 = p.intersections(e3);

    REQUIRE(i3.size() == 2);
    CHECK(std::min(i3[0].x, i3[1].x) == 1.0_e13);
    CHECK(std::max(i3[0].x, i3[1].x) == 2.0_e13);

    CHECK(std::min(i3[0].y, i3[1].y) == 1.5_e13);
    CHECK(std::max(i3[0].y, i3[1].y) == 1.5_e13);
  }
}

TEST_CASE("geometry/2d/barrel/polygon/intersection/math") {
  std::ifstream in("polygon.test");

  if (!in) {
    WARN("cannot open file `polygon.test'");
    return;
  }

  int n_events;
  in >> n_events;

  Polygon<4> poly;

  for (int i = 0; i < 4; ++i) {
    double x, y;
    in >> x >> y;
    Polygon<4>::Point p(x, y);
    poly.push_back(p);
  }

  for (int i = 0; i < n_events; i++) {
    double x, y, phi;
    in >> x >> y >> phi;

    double a, b, c;
    in >> a >> b >> c;

    size_t n_iters;
    in >> n_iters;

    Polygon<4>::Event event(x, y, phi);
    bool intersects = n_iters > 0;

    CHECKED_IF(poly.intersects(event) == intersects) {
      auto inters = poly.intersections(event);

      CHECKED_IF(inters.size() == n_iters) {

        if (n_iters > 0) {
          double ix, iy;
          Polygon<4>::Intersections m_inters;

          for (size_t j = 0; j < n_iters; ++j) {
            in >> ix >> iy;
            Polygon<4>::Point p(ix, iy);
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
