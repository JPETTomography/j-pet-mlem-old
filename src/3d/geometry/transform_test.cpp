#include "util/test.h"

#include "3d/geometry/ray.h"
#include "3d/geometry/point.h"

#include "3d/geometry/transform.h"

#include "common/types.h"

using Vector = PET3D::Vector<F>;
using Point = PET3D::Point<F>;

TEST_CASE("3d/transform") {

  SECTION("vector/rotation/1") {
    Vector src(0.1f, 0.2f, 0.3f);
    auto dest = rotate(src, F(M_PI / 3.0), Vector(0, 0, 1.0));

    CHECK(dest == VApprox(Vector(-0.123205080757, 0.186602540378, 0.3)));
  }

  SECTION("Vector/rotation/2") {
    Vector src(0.1f, 0.2f, 0.3f);
    auto dest = rotate(src, F(M_PI / 3.0), Vector(.2, .4, 1).normalized());

    CHECK(dest ==
          VApprox(Vector(0.02008778013, 0.198289443268, 0.316666666667)));
  }

  SECTION("point/rotation/1") {
    Point src(0.1f, 0.2f, 0.3f);
    auto dest = rotate(src, F(M_PI / 3.0), Vector(0, 0, 1.0));

    CHECK(dest == VApprox(Point(-0.123205080757, 0.186602540378, 0.3)));
  }

  SECTION("point/rotation/2") {
    Point src(0.1f, 0.2f, 0.3f);
    auto dest = rotate(src, F(M_PI / 3.0), Vector(.2, .4, 1).normalized());

    CHECK(dest ==
          VApprox(Point(0.02008778013, 0.198289443268, 0.316666666667)));
  }

  SECTION("point/rotation/2") {
    Point src(0.1f, 0.2f, 0.3f);
    auto dest = rotate(
        src, F(M_PI / 3.0), Vector(.2, .4, 1).normalized(), Point(0.2, .4, .7));

    CHECK(dest == VApprox(Point(0.139956, 0.200855, 0.291667)));
  }
}
