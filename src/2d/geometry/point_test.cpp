#include "util/test.h"

#include "point.h"
#include "transformation.h"
#include "common/types.h"

TEST("Transforming and comparing points") {

  using Point = PET2D::Point<F>;
  using Transformation = PET2D::Transformation<F>;
  using Vector = typename Point::Vector;

  Point p(0, 1);
  Point q(0, 1);

  REQUIRE(p == q);
  REQUIRE(p.approx_equal(q));

  auto t1 = Transformation(M_PI / 3);
  auto t2 = Transformation(M_PI / 3);
  auto t3 = Transformation(2 * M_PI / 3);

  REQUIRE(t2(t1(p)).approx_equal(t3(p)));
}
