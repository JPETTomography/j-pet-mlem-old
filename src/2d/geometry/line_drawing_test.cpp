#include <iostream>

#include "util/test.h"

#include "line_drawing.h"

TEST("2d/geometry/line_drawing") {
  PET2D::PixelGrid<float, int> grid(
      128, 128, 0.005, PET2D::Point<float>(-64 * 0.005, -64 * 0.005));

  {
    PET2D::Point<float> start(0.001, 0.001);
    PET2D::Point<float> end(0.007, 0.003);
    using Container = std::vector<PET2D::Pixel<int>>;
    Container pixels;
    PET2D::draw_line(start, end, grid, std::back_inserter(pixels));
#if THIS_IS_NOT_A_TEST
    // FIXME: this is NOT a test
    std::cout << "----\n";
    for (PET2D::Pixel<int> p : pixels) {
      std::cout << p.x << " " << p.y << "\n";
    }
#endif
  }

  {
    PET2D::Point<float> start(0.001, 0.001);
    PET2D::Point<float> end(0.001, -0.010);
    using Container = std::vector<PET2D::Pixel<int>>;
    Container pixels;
    PET2D::draw_line(start, end, grid, std::back_inserter(pixels));
#if THIS_IS_NOT_A_TEST
    // FIXME: this is NOT a test
    std::cout << "----\n";
    for (PET2D::Pixel<int> p : pixels) {
      std::cout << p.x << " " << p.y << "\n";
    }
#endif
  }

  {
    PET2D::Point<float> start(0.001, 0.001);
    PET2D::Point<float> end(0.020, 0.001);
    using Container = std::vector<PET2D::Pixel<int>>;
    Container pixels;
    PET2D::draw_line(start, end, grid, std::back_inserter(pixels));
#if THIS_IS_NOT_A_TEST
    // FIXME: this is NOT a test
    std::cout << "----\n";
    for (PET2D::Pixel<int> p : pixels) {
      std::cout << p.x << " " << p.y << "\n";
    }
#endif
  }
}
