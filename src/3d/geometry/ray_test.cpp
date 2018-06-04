#include "util/test.h"

#include "3d/geometry/ray.h"
#include "3d/geometry/point.h"

#include "common/types.h"

using Vector = PET3D::Vector<F>;
using Point = PET3D::Point<F>;

TEST_CASE("3d/geometry/ray") {
  using namespace ray_tracing;
  Point p(F(0.1), F(.2), F(0.3));
  Vector d(F(.3), F(.1), F(.2));

  Ray<F> ray(p, d);

  SECTION("init") {
    CHECK(ray.p == VApprox(p));
    CHECK(ray.d == VApprox(d));
  }
}

TEST_CASE("3d/geometry/box") {
  using namespace ray_tracing;

  SECTION("init/AAB") {
    Box<F> box = Box<F>::AAB(Point(0, 0, 0), Point(1, 1, 1));

    CHECK(box.center == VApprox(Point(0.5, 0.5, 0.5)));

    CHECK(box.a_u == VApprox(Vector::e_x()));
    CHECK(box.a_v == VApprox(Vector::e_y()));
    CHECK(box.a_w == VApprox(Vector::e_z()));

    CHECK(box.h_u == Approx(0.5));
    CHECK(box.h_v == Approx(0.5));
    CHECK(box.h_w == Approx(0.5));
  }
}

TEST_CASE("3d/geometry/ray_box") {
  using namespace ray_tracing;

  SECTION("intersect 1") {
    Box<F> box = Box<F>::AAB(Point(0, 0, 0), Point(1, 1, 1));
    Ray<F> ray(Point(0.5, 0.5, 0.5), Vector(0, 1, 0).normalized());

    auto intersection = intersect(ray, box);

    CHECK(intersection.intersected);
    F t;
    if (intersection.t_min > 0)
      t = intersection.t_min;
    else
      t = intersection.t_max;
    CHECK(t == Approx(0.5));

    CHECK(ray(t) == VApprox(Point(0.5, 1.0, 0.5)));
  }

  SECTION("intersect 2") {
    Box<F> box = Box<F>::AAB(Point(0, 0, 0), Point(1, 1, 1));
    Ray<F> ray(Point(0.5, 0.5, 0.5), Vector(0, 1, 0.1).normalized());

    auto intersection = intersect(ray, box);

    CHECK(intersection.intersected);
    F t;
    if (intersection.t_min > 0)
      t = intersection.t_min;
    else
      t = intersection.t_max;
    CHECK(t == Approx(0.50249378));

    CHECK(ray(t) == VApprox(Point(0.5, 1.0, 0.55)));
  }

  SECTION("intersect/rotated") {
    Box<F> box = Box<F>::rotate(Box<F>::AAB(Point(0, 0, 0), Point(1, 1, 1)),
                                F(M_PI / 3),
                                Vector(0.5, .75, 1).normalized(),
                                Point(0, 0, 0));

    Ray<F> ray(Point(0.5, 0.5, 0.5), Vector(0, 1, 0.1).normalized());

    auto intersection = intersect(ray, box);

    CHECK(intersection.intersected);
    F t;
    if (intersection.t_min > 0)
      t = intersection.t_min;
    else
      t = intersection.t_max;
    CHECK(t == Approx(0.725833346158));

    CHECK(ray(t) == VApprox(Point(0.5, 1.22223117332, 0.572223117332)));
  }
}
