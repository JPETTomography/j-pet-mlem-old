#include "util/test.h"

#include "transformation.h"

#include "vector.h"
#include "point.h"

#include "common/types.h"

TEST("2d transformation") {
  using Vector = PET2D::Vector<F>;
  using Point = PET2D::Point<F>;
  using Transformation = PET2D::Transformation<F>;

  SECTION("Identity translation") {
    Point p(1, 0.5);
    Transformation id;
    Point pt = id(p);
    CHECK(pt.x == Approx(1.0));
    CHECK(pt.y == Approx(0.5));
  }

  SECTION("Single translation") {
    Point p(1, 0.5);
    Transformation translate(Vector(0.3, 0.5));
    Point pt = translate(p);
    CHECK(pt.x == Approx(1.3));
    CHECK(pt.y == Approx(1.0));

    Vector v(1.0, 0.5);
    auto vt = translate(v);
    CHECK(vt.x == Approx(1));
    CHECK(vt.y == Approx(0.5));
  }

  SECTION("Single rotation and translation") {
    Point p(1, 0.5);
    Transformation transform(M_PI / 6, Vector(0.3, 0.5));
    Point pt = transform(p);
    CHECK(pt.x == Approx(0.916025));
    CHECK(pt.y == Approx(1.43301));

    Vector v(1.0, 0.5);
    auto vt = transform(v);
    CHECK(vt.x == Approx(0.616025));
    CHECK(vt.y == Approx(0.933013));
  }

  SECTION("two transform composition") {
    Point p(1, 0.5);
    Transformation transform1(M_PI / 6, Vector(0.3, 0.5));
    Transformation transform2(M_PI / 4, Vector(-0.5, 1.0));
    auto transform = transform2 * transform1;
    Point pt = transform(p);
    CHECK(pt.x == Approx(-0.865565));
    CHECK(pt.y == Approx(2.66102));

    Vector v(1.0, 0.5);
    auto vt = transform(v);
    CHECK(vt.x == Approx(-0.224144));
    CHECK(vt.y == Approx(1.09534));
  }

  SECTION("copy constructor") {
    Transformation transform1(M_PI / 6, Vector(0.3, 0.5));
    Transformation transform2 = transform1;
    CHECK(transform1.rotation == transform2.rotation);
    CHECK(transform1.translation == transform2.translation);
  }
}
