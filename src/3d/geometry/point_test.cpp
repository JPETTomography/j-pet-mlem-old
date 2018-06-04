#include "util/test.h"

#include "3d/geometry/point.h"

#include "common/types.h"

using Point = PET3D::Point<F>;
using Vector = PET3D::Vector<F>;

TEST("3d/geometry/point/init", "point construction") {
  Point p(1, 2, 3);

  CHECK(p.x == 1.0_e7);
  CHECK(p.y == 2.0_e7);
  CHECK(p.z == 3.0_e7);
}

TEST("3d/geometry/point/arithmetic assignemt", "point arithmetic assignment") {
  {
    Point p(1, 2, 3);
    Vector v(0.1, 0.2, 0.3);
    p += v;
    CHECK(p.x == 1.1_e7);
    CHECK(p.y == 2.2_e7);
    CHECK(p.z == 3.3_e7);
  }
  {
    Point p(1, 2, 3);
    Vector v(0.1, 0.2, 0.3);
    p -= v;
    CHECK(p.x == 0.9_e7);
    CHECK(p.y == 1.8_e7);
    CHECK(p.z == 2.7_e7);
  }
}

TEST("3d/geometry/point/difference", "point differencet") {
  Point p1(1, 2, 3);
  Point p2(0.1, 0.2, 0.3);

  Vector v = p1 - p2;

  CHECK(v.x == 0.9_e7);
  CHECK(v.y == 1.8_e7);
  CHECK(v.z == 2.7_e7);
}

TEST("3d/geometry/point/distance", "point distances") {
  Point p1(1, 2, 3);

  CHECK(p1.distance_from_origin() == Approx(std::sqrt(14)).epsilon(1e-7));
  CHECK(p1.distance_from_origin2() == Approx(14).epsilon(1e-7));
}

TEST("3d/geometry/point/nearest_distance", "point nearest distance") {
  Point p1(1, 2, 3);
  Point p2(1, 3, 3.3);
  Point p3(1, 2.2, 4);

  CHECK(p1.nearest_distance(p2, p3) == Approx(std::sqrt(1.04)).epsilon(1e-7));
}
