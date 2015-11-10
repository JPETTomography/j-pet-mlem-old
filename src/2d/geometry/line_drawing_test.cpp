#include <iostream>

#include "util/test.h"

#include "line_drawing.h"

#include "common/types.h"

using Grid = PET2D::PixelGrid<F, S>;
using Point = PET2D::Point<F>;
using Pixel = PET2D::Pixel<S>;

TEST("2d/geometry/line_drawing") {
  Grid grid(128, 128, 0.005, Point(-64 * 0.005, -64 * 0.005));
  {
    Point start(0.001, 0.001);
    Point end(0.007, 0.003);
    using Container = std::vector<Pixel>;
    Container pixels;
    PET2D::draw_line(start, end, grid, std::back_inserter(pixels));
#if THIS_IS_NOT_A_TEST
    // FIXME: this is NOT a test
    std::cout << "----\n";
    for (Pixel p : pixels) {
      std::cout << p.x << ' ' << p.y << "\n";
    }
#endif
  }
  {
    Point start(0.001, 0.001);
    Point end(0.001, -0.010);
    using Container = std::vector<Pixel>;
    Container pixels;
    PET2D::draw_line(start, end, grid, std::back_inserter(pixels));
#if THIS_IS_NOT_A_TEST
    // FIXME: this is NOT a test
    std::cout << "----\n";
    for (Pixel p : pixels) {
      std::cout << p.x << ' ' << p.y << "\n";
    }
#endif
  }
  {
    Point start(0.001, 0.001);
    Point end(0.020, 0.001);
    using Container = std::vector<Pixel>;
    Container pixels;
    PET2D::draw_line(start, end, grid, std::back_inserter(pixels));
#if THIS_IS_NOT_A_TEST
    // FIXME: this is NOT a test
    std::cout << "----\n";
    for (Pixel p : pixels) {
      std::cout << p.x << ' ' << p.y << "\n";
    }
#endif
  }
}
