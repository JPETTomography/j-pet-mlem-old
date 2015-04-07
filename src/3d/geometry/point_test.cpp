#include "util/test.h"

#include "3d/geometry/point.h"

using namespace PET3D;

TEST("3d/geometry/point/init", "point construction") {
  using Point = PET3D::Point<float>;
  using Vector = PET3D::Vector<float>;

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
    p+=v;
    CHECK(p.x == 1.1_e7);
    CHECK(p.y == 2.2_e7);
    CHECK(p.z == 3.3_e7);
  }

    {
      Point p(1.0f, 2.0f, 3.0f);
      Vector v(0.1f, 0.2f, 0.3f);
      p-=v;
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

    Vector v=p1-p2;

    CHECK(v.x == 0.9_e7);
    CHECK(v.y == 1.8_e7);
    CHECK(v.z == 2.7_e7);

}
