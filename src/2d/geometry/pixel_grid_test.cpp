#include "util/test.h"

#include "2d/geometry/pixel_grid.h"

#include "common/types.h"

using Point = PET2D::Point<F>;
using Pixel = PET2D::Pixel<S>;
using PixelGrid = PET2D::PixelGrid<F, S>;

TEST("2d/geometry/pixel_grid/pixel_at") {
  PixelGrid grid(128, 64, 0.005, Point(-0.005 * 64, -0.005 * 32));
  {
    Pixel p = grid.pixel_at(Point(0, 0));
    REQUIRE(p.x == 64);
    REQUIRE(p.y == 32);
  }
  {
    Pixel p = grid.pixel_at(Point(-0.006, 0.019));
    REQUIRE(p.x == 62);
    REQUIRE(p.y == 35);
  }
}

TEST("2d/geometry/pixel_grid/point_at") {
  PixelGrid grid(128, 64, 0.005, Point(-0.005 * 64, -0.005 * 32));
  {
    Point p = grid.lower_left_at(64, 32);
    REQUIRE(p.x == 0.0_e7);
    REQUIRE(p.y == 0.0_e7);
  }
  {
    Point p = grid.center_at(64, 32);
    REQUIRE(p.x == 0.0025_e7);
    REQUIRE(p.y == 0.0025_e7);
  }
  {
    Point p = grid.lower_left_at(66, 35);
    REQUIRE(p.x == 0.010_e7);
    REQUIRE(p.y == 0.015_e7);
  }
  {
    Point p = grid.center_at(66, 35);
    REQUIRE(p.x == 0.0125_e7);
    REQUIRE(p.y == 0.0175_e7);
  }
}
