#include "catch.hpp"

#include <cmath>

#include "pixel_grid.h"

TEST_CASE("pixel_grid", "pixel grid") {

  Point<> ll(-250, -200);
  Point<> ur(100, 50);

  int nx = 70;
  int ny = 125;

  PixelGrid<> grid(ll, ur, nx, ny);

  double dx = 5.0;
  double dy = 2.0;

  SECTION("point_index_set", "") {

    double x = 1.0;
    double y = 2.8;
    Point<> p(x, y);

    REQUIRE(x == p.x);
    REQUIRE(y == p.y);

    int ix = 10;
    int iy = 28;
    Index<> ind(ix, iy);

    REQUIRE(ix == ind.x);
    REQUIRE(iy == ind.y);
  }

  SECTION("set", "") {

    REQUIRE(nx == grid.nx());
    REQUIRE(ny == grid.ny());

    REQUIRE(dx == grid.dx());
    REQUIRE(dy == grid.dy());
  }

  SECTION("index", "") { REQUIRE((23 * 70 + 13) == grid.index(13, 23)); }

  SECTION("center", "") {
    {
      Point<> c = grid.center(0, 0);

      REQUIRE(-247.5 == c.x);
      REQUIRE(-199.0 == c.y);
    }
    {
      Point<> c = grid.center(50, 70);

      REQUIRE(2.5 == c.x);
      REQUIRE(-59.0 == c.y);
    }
  }

  SECTION("in", "") {
    {
      Point<> c = grid.center(0, 0);

      REQUIRE(-247.5 == c.x);
      REQUIRE(-199.0 == c.y);

      Index<> ind = grid.in(c);
      REQUIRE(0 == ind.x);
      REQUIRE(0 == ind.y);
    }
    {
      Point<> c = grid.center(50, 70);

      REQUIRE(2.5 == c.x);
      REQUIRE(-59.0 == c.y);

      Index<> ind = grid.in(c);
      REQUIRE(50 == ind.x);
      REQUIRE(70 == ind.y);
    }
    {
      Point<> c = grid.center(0, 0);

      REQUIRE(-247.5 == c.x);
      REQUIRE(-199.0 == c.y);

      Index<> ind = grid.in(Point<>(c.x + 1.75, c.y - 0.999));
      REQUIRE(0 == ind.x);
      REQUIRE(0 == ind.y);

    }
    {
      Point<> c = grid.center(50, 70);

      REQUIRE(2.5 == c.x);
      REQUIRE(-59.0 == c.y);
      Index<> ind = grid.in(Point<>(c.x - 2.4999, c.y + 0.999));
      REQUIRE(50 == ind.x);
      REQUIRE(70 == ind.y);
    }
    {
      Point<> c = grid.center(50, 70);

      REQUIRE(2.5 == c.x);
      REQUIRE(-59.0 == c.y);
      Index<> ind = grid.in(Point<>(c.x - 2.4999, c.y + 1.0001));
      REQUIRE(50 == ind.x);
      REQUIRE(70 != ind.y);
    }
  }

  SECTION("add", "") {

    grid.add(0, 0, 0.5);

    REQUIRE(0.5 == grid[0]);
    grid.add(0, 0, 0.76);
    REQUIRE(1.26 == grid[0]);

    grid.add(13, 23);
    REQUIRE(1.0 == grid(23 * 70 + 13));
  }

  SECTION("insert", "") {
    {
      Point<> c = grid.center(0, 0);

      REQUIRE(-247.5 == c.x);
      REQUIRE(-199.0 == c.y);

      grid.insert(-247.5, -199.0, 0.5);

      REQUIRE(0.5 == grid(0, 0));
    }
    {
      Point<> c = grid.center(50, 70);

      REQUIRE(2.5 == c.x);
      REQUIRE(-59.0 == c.y);

      grid.insert(2.5, -59.0, 0.5);

      REQUIRE(0.5 == grid(50, 70));
    }
  }
}
