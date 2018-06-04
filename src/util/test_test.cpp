#include "util/test.h"

#include "3d/geometry/vector.h"
#include "3d/geometry/point.h"

#include "common/types.h"

using Vector = PET3D::Vector<F>;
using Point = PET3D::Point<F>;

TEST_CASE("VApprox") {

  //  CHECK(Vector(0, 0, 0) == VApprox(Vector(0, 0, 1)));
  //  CHECK(Point(0, 0, 0) == VApprox(Point(0, 0, 1)));

  //  CHECK(VApprox(Vector(0, 0, 1)) == Vector(0, 0, 0));
  //  CHECK(VApprox(Point(0, 0, 1)) == Point(0, 0, 0));

  CHECK(Vector(0, 0, 1) == VApprox(Vector(0, 0, 1)));
  CHECK(Point(0, 0, 1) == VApprox(Point(0, 0, 1)));
}
