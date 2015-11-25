#include <iostream>
#include <fstream>

#include "util/test.h"

#include "polygon.h"

#include "common/types.h"

using Point = PET2D::Point<F>;
using Polygon = PET2D::Polygon<4, F>;
using Polygon3 = PET2D::Polygon<3, F>;

TEST("2d/geometry/polygon/json") {
  Polygon ps;
  ps.emplace_back(1, 1);
  ps.emplace_back(2, 1);
  ps.emplace_back(2, 2);
  ps.emplace_back(1, 2);

  json j(ps);
  REQUIRE(j.dump() == "[[1,1],[2,1],[2,2],[1,2]]");
}

TEST("2d/geometry/polygon/center") {
  Polygon ps;
  ps.emplace_back(1, 1);
  ps.emplace_back(2, 1);
  ps.emplace_back(2, 2);
  ps.emplace_back(1, 2);

  CHECK(ps.center() == Point(1.5, 1.5));

  Polygon3 pt;
  pt.emplace_back(1, 1);
  pt.emplace_back(2, 1);
  pt.emplace_back(1.5, 2);

  CHECK(pt.center() == Point(1.5, 4. / 3.));
}

TEST("2d/geometry/polygon/intersection") {
  Polygon p;
  p.emplace_back(1, 1);
  p.emplace_back(2, 1);
  p.emplace_back(2, 2);
  p.emplace_back(1, 2);

  Polygon::Event e1(0, 0.5, M_PI_4);
  Polygon::Event e2(0, 1.5, M_PI_4);
  Polygon::Event e3(0, 1.5, 0);

  CHECK(true == p.intersects(e1));
  CHECK(false == p.intersects(e2));
  CHECK(true == p.intersects(e3));

  auto i1 = p.intersections(e1);

  REQUIRE(i1.size() == 2);
  CHECK(std::min(i1[0].x, i1[1].x) == 1._e13);
  CHECK(std::max(i1[0].x, i1[1].x) == 1.5_e13);

  CHECK(std::min(i1[0].y, i1[1].y) == 1.5_e13);
  CHECK(std::max(i1[0].y, i1[1].y) == 2._e13);

  auto i3 = p.intersections(e3);

  REQUIRE(i3.size() == 2);
  CHECK(std::min(i3[0].x, i3[1].x) == 1._e13);
  CHECK(std::max(i3[0].x, i3[1].x) == 2._e13);

  CHECK(std::min(i3[0].y, i3[1].y) == 1.5_e13);
  CHECK(std::max(i3[0].y, i3[1].y) == 1.5_e13);
}

TEST("2d/geometry/polygon/intersection/math") {
  std::ifstream in("test_input/polygon_test.tab");

  if (!in) {
    WARN(
        "cannot open file `test_input/polygon_test.tab', "
        "evaluate `math/polygon_test.nb'");
    return;
  }

  int n_events;
  in >> n_events;

  Polygon poly;

  for (int i = 0; i < 4; ++i) {
    double x, y;
    in >> x >> y;
    Point p(x, y);
    poly.push_back(p);
  }

  for (int i = 0; i < n_events; i++) {
    double x, y, phi;
    in >> x >> y >> phi;

    double a, b, c;
    in >> a >> b >> c;

    size_t n_iters;
    in >> n_iters;

    Polygon::Event event(x, y, phi);
    bool intersects = n_iters > 0;

    CHECKED_IF(poly.intersects(event) == intersects) {
      auto inters = poly.intersections(event);

      CHECKED_IF(inters.size() == n_iters) {

        if (n_iters > 0) {
          double in_x, in_y;
          Polygon::Intersections m_inters;

          for (size_t j = 0; j < n_iters; ++j) {
            in >> in_x >> in_y;
            Point p(in_x, in_y);
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

TEST("2d/geometry/polygon/is_inside") {
  Polygon ps;

  ps.emplace_back(0, 0);
  ps.emplace_back(0, 1);
  ps.emplace_back(1, 1);
  ps.emplace_back(1, 0);

  CHECK(ps.contains(Point(0.5, 0.5)) == true);
  CHECK(ps.contains(Point(0.7, 1.1)) == false);
  CHECK(ps.contains(Point(1.0e-12, 0.5)) == true);
  CHECK(ps.contains(Point(-1.0e-12, 0.5)) == false);
  CHECK(ps.contains(Point(-.01, .6)) == false);
}
