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
