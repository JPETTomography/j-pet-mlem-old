#include "util/test.h"

#include "3d/geometry/point.h"

using namespace PET3D;

TEST("3d/geometry/point/init", "point construction") {
  using Point = PET3D::Point<float>;


  Point p(1.0f, 2.0f, 3.0f);

  CHECK(p.x == 1.0_e7);
  CHECK(p.y == 2.0_e7);
  CHECK(p.z == 3.0_e7);
}

TEST("3d/geometry/point/arithmetic assignemt", "point arithmetic assignment") {
  using Point = PET3D::Point<float>;
  using Vector = PET3D::Vector<float>;

  {
    Point p(1.0f, 2.0f, 3.0f);
    Vector v(0.1f, 0.2f, 0.3f);
    p += v;
    CHECK(p.x == 1.1_e7);
    CHECK(p.y == 2.2_e7);
    CHECK(p.z == 3.3_e7);
  }

  {
    Point p(1.0f, 2.0f, 3.0f);
    Vector v(0.1f, 0.2f, 0.3f);
    p -= v;
    CHECK(p.x == 0.9_e7);
    CHECK(p.y == 1.8_e7);
    CHECK(p.z == 2.7_e7);
  }
}

TEST("3d/geometry/point/difference", "point differencet") {
  using Point = PET3D::Point<float>;
  using Vector = PET3D::Vector<float>;

  Point p1(1.0f, 2.0f, 3.0f);
  Point p2(0.1f, 0.2f, 0.3f);

  Vector v = p1 - p2;

  CHECK(v.x == 0.9_e7);
  CHECK(v.y == 1.8_e7);
  CHECK(v.z == 2.7_e7);
}

TEST("3d/geometry/point/distance", "point distances") {
  using Point = PET3D::Point<float>;

  Point p1(1.0f, 2.0f, 3.0f);

  CHECK(p1.distance_from_origin() == Approx(std::sqrt(14.0f)).epsilon(1e-7));
  CHECK(p1.distance_from_origin2() == Approx(14.0f).epsilon(1e-7));
}

TEST("3d/geometry/point/nearest_distance", "point nearest distance") {
  using Point = PET3D::Point<float>;

  Point p1(1.0f, 2.0f, 3.0f);
  Point p2(1.0f, 3.0f, 3.30f);
  Point p3(1.0f, 2.2f, 4.0f);

  CHECK(p1.nearest_distance(p2, p3) == Approx(std::sqrt(1.04f)).epsilon(1e-7));
}
