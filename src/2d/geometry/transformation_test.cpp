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
    Transformation tranform(M_PI / 6, Vector(0.3, 0.5));
    Point pt = tranform(p);
    CHECK(pt.x == Approx(0.916025));
    CHECK(pt.y == Approx(1.43301));
  }

  SECTION("two transform composition") {
    Point p(1, 0.5);
    Transformation transform1(M_PI / 6, Vector(0.3, 0.5));
    Transformation transform2(M_PI / 4, Vector(-0.5, 1.0));
    auto transform = transform2 * transform1;
    Point pt = transform(p);
    CHECK(pt.x == Approx(-0.865565));
    CHECK(pt.y == Approx(2.66102));
  }
}
