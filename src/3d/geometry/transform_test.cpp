#include "util/test.h"

#include "3d/geometry/ray.h"
#include "3d/geometry/point.h"

#include "3d/geometry/transform.h"

#include "common/types.h"

using Vector = PET3D::Vector<F>;
using Point = PET3D::Point<F>;

TEST_CASE("3d/transform") {
  Vector src(0.1f, 0.2f, 0.3f);
  auto dest = rotate(src, F(M_PI / 3.0), Vector(0, 0, 1.0));

  CHECK(dest == VApprox(Vector(0.326794919243, 0.966025403784, 0.3)));
}
