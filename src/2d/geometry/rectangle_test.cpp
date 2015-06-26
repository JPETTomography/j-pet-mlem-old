#include <fstream>

#include "util/test.h"

#include "rectangle.h"

TEST("2d/geometry/rectangle") {
  using Point = PET2D::Point<float>;

  PET2D::Rectangle<float> r(1, 2, 3, 4);

  REQUIRE(r.area == 48.0_e7);
  REQUIRE(r.contains(Point(1, 2)));
  REQUIRE(r.contains(Point(1, 5.9)));
  REQUIRE(r.contains(Point(-1.99, 2)));

  REQUIRE_FALSE(r.contains(Point(1, 6.1)));
  REQUIRE_FALSE(r.contains(Point(-2.1, 2)));
}

TEST("2d/geometry/rectangle/point_generator") {
  using Point = PET2D::Point<float>;

  PET2D::Rectangle<float> r(1, 2, 3, 4);

  PET2D::RectanglePointGenerator<float> gen(r);
  std::mt19937_64 rng;

  std::ofstream out("test_output/random_rectangle_point.txt");
  for (int i = 0; i < 100; i++) {
    auto p = gen(rng);
    REQUIRE(r.contains(p));
    out << p << "\n";
  }
}
