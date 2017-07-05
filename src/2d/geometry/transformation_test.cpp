#include "util/test.h"

#include "transformation.h"

#include "vector.h"
#include "point.h"

#include "common/types.h"

TEST("2d transformation") {
  using Vector = PET2D::Vector<F>;
  using Point = PET2D::Point<F>;
  using Transformation = PET2D::Transformation<F>;

  SECTION("Single translation") {
    Point p(1, 0.5);
    Transformation translate(Vector(0.3, 0.5));
    Point pt = translate(p);
    CHECK(pt.x == Approx(1.3));
    CHECK(pt.y == Approx(1.0));
  }

  SECTION("Single rotation and translation") {
    Point p(1, 0.5);
    Transformation tranaform(M_PI / 6, Vector(0.3, 0.5));
    Point pt = tranaform(p);
    CHECK(pt.x == Approx(0.916025));
    CHECK(pt.y == Approx(1.43301));
  }
}
