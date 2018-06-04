#include "util/test.h"

#include <cmath>

#include "ellipse.h"

#include "common/types.h"

using Ellipse = PET2D::Ellipse<F>;

TEST("2d/geometry/ellipse") {}

TEST("2d/geometry/ellipse/write_read") {
  Ellipse ellipse(1, 2, 3, 4, 90);
  {
    auto fn = "ellipse_test.txt"_temp;
    std::ofstream out(fn);
    out << ellipse;
    out.close();

    std::ifstream in(fn);
    REQUIRE(in);
    Ellipse ellipse_read(in);
    in.close();

    CHECK(ellipse.center == ellipse_read.center);
    CHECK(ellipse.a == ellipse_read.a);
    CHECK(ellipse.b == ellipse_read.b);
    CHECK(ellipse.angle == ellipse_read.angle);

    REQUIRE(std::remove(fn.c_str()) == 0);
  }
}
