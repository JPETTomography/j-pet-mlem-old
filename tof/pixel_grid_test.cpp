#include <cmath>
#include <gtest/gtest.h>
#include "pixel_grid.h"

class pixel_grid_test : public ::testing::Test {
protected:

  pixel_grid_test():
    ll(Point<>(-250, -200)), ur(Point<>(100, 50)), nx(70), ny(125), grid(ll, ur, nx, ny) {};
  virtual void SetUp() {
    dx=5;
    dy=2;
  }

#if 1
  int nx;
  int ny;

  Point<> ll;
  Point<> ur;
  PixelGrid<> grid;

  double dx, dy;
#endif
};
TEST_F(pixel_grid_test, point_index_set_test) {
  double x=1.0;
  double y=2.8;
  Point<> p(x, y);

  ASSERT_DOUBLE_EQ(x, p.x);
  ASSERT_DOUBLE_EQ(y, p.y);

  int ix=10;
  int iy=28;
  Index<> ind(ix, iy);

  ASSERT_EQ(ix, ind.x);
  ASSERT_EQ(iy, ind.y);
}

TEST_F(pixel_grid_test, set_test) {

  ASSERT_EQ(nx, grid.nx());
  ASSERT_EQ(ny, grid.ny());

  ASSERT_DOUBLE_EQ(dx, grid.dx());
  ASSERT_DOUBLE_EQ(dy, grid.dy());
}
TEST_F(pixel_grid_test, index_test) {
  ASSERT_EQ(23*70+13, grid.index(13, 23));
}

TEST_F(pixel_grid_test, center_test) {
  {
    Point<> c=grid.center(0, 0);

    ASSERT_DOUBLE_EQ(-247.5, c.x);
    ASSERT_DOUBLE_EQ(-199.0, c.y);
  }
  {
    Point<> c=grid.center(50, 70);

    ASSERT_DOUBLE_EQ(  2.5, c.x);
    ASSERT_DOUBLE_EQ(-59.0, c.y);
  }
}
TEST_F(pixel_grid_test, in_test) {
  {
    Point<> c=grid.center(0, 0);

    ASSERT_DOUBLE_EQ(-247.5, c.x);
    ASSERT_DOUBLE_EQ(-199.0, c.y);

    Index<> ind=grid.in(c);
    ASSERT_EQ(0, ind.x);
    ASSERT_EQ(0, ind.y);

  }
  {
    Point<> c=grid.center(50, 70);

    ASSERT_DOUBLE_EQ(  2.5, c.x);
    ASSERT_DOUBLE_EQ(-59.0, c.y);

    Index<> ind=grid.in(c);
    ASSERT_EQ(50, ind.x);
    ASSERT_EQ(70, ind.y);

  }

    {
    Point<> c=grid.center(0, 0);

    ASSERT_DOUBLE_EQ(-247.5, c.x);
    ASSERT_DOUBLE_EQ(-199.0, c.y);

    Index<> ind=grid.in(Point<>(c.x+1.75, c.y-0.999) );
    ASSERT_EQ(0, ind.x);
    ASSERT_EQ(0, ind.y);

  }
  {
    Point<> c=grid.center(50, 70);

    ASSERT_DOUBLE_EQ(  2.5, c.x);
    ASSERT_DOUBLE_EQ(-59.0, c.y);
    Index<> ind=grid.in(Point<>(c.x-2.4999, c.y+0.999) );
    ASSERT_EQ(50, ind.x);
    ASSERT_EQ(70, ind.y);

  }

    {
    Point<> c=grid.center(50, 70);

    ASSERT_DOUBLE_EQ(  2.5, c.x);
    ASSERT_DOUBLE_EQ(-59.0, c.y);
    Index<> ind=grid.in(Point<>(c.x-2.4999, c.y+1.0001) );
    ASSERT_EQ(50, ind.x);
    ASSERT_NE(70, ind.y);

  }
}
TEST_F(pixel_grid_test, add_test) {
  grid.add(0, 0, 0.5);

  ASSERT_DOUBLE_EQ(0.5, grid[0]);
  grid.add(0, 0, 0.76);
  ASSERT_DOUBLE_EQ(1.26, grid[0]);

  grid.add(13, 23);
  ASSERT_EQ(1.0, grid(23*70+13));
}
TEST_F(pixel_grid_test, insert_test) {
  {
    Point<> c=grid.center(0, 0);

    ASSERT_DOUBLE_EQ(-247.5, c.x);
    ASSERT_DOUBLE_EQ(-199.0, c.y);

    grid.insert(-247.5, -199.0, 0.5);

    ASSERT_DOUBLE_EQ(0.5, grid(0, 0));
  }
  {
    Point<> c=grid.center(50, 70);

    ASSERT_DOUBLE_EQ(  2.5, c.x);
    ASSERT_DOUBLE_EQ(-59.0, c.y);

    grid.insert(2.5, -59.0, 0.5);

    ASSERT_DOUBLE_EQ(0.5, grid(50, 70));
  }
}
