#include <fstream>

#include "util/test.h"

#include "mathematica_graphics.h"

#include "2d/barrel/square_detector.h"
#include "2d/barrel/barrel_builder.h"
#include "2d/geometry/line_segment.h"
#include "2d/geometry/pixel_grid.h"

using F = float;
using S = int;

using Detector = PET2D::Barrel::SquareDetector<F>;
using BarrelBuilder = PET2D::Barrel::BarrelBuilder<Detector, S>;
using Scanner = BarrelBuilder::BigBarrel;
using MathematicaGraphics = Common::MathematicaGraphics<F>;

TEST("common/mathematica_graphics/detector") {
  std::ofstream out("test_output/graph_detector.m");
  if (!out) {
    FAIL("cannot open file");
  }

  MathematicaGraphics graphics(out);

  Detector detector(0.007, 0.019, 0);

  graphics.add(detector);
}

TEST("common/mathematica_graphics/big_barrel") {

  auto scanner = BarrelBuilder::make_big_barrel();

  std::ofstream out("test_output/graph_scanner.m");
  if (!out) {
    FAIL("cannot open file");
  }

  MathematicaGraphics graphics(out);

  graphics.add(scanner);
}

TEST("common/mathematica_graphics/big_barrel/lor") {

  auto scanner = BarrelBuilder::make_big_barrel();

  std::ofstream out("test_output/graph_lor.m");
  if (!out) {
    FAIL("cannot open file");
  }

  MathematicaGraphics graphics(out);

  graphics.add(scanner);
  graphics.add(scanner, PET2D::Barrel::LOR<int>(65, 0));
}

TEST("common/mathematica_graphics/big_barrel/segment") {

  auto scanner = BarrelBuilder::make_big_barrel();
  using Point = PET2D::Point<F>;

  std::ofstream out("test_output/graph_segment.m");
  if (!out) {
    FAIL("cannot open file");
  }

  MathematicaGraphics graphics(out);

  graphics.add(scanner);

  PET2D::LineSegment<F> segment(Point(-0.400, 0), Point(0, 0.400));

  graphics.add(segment);
}

TEST("common/mathematica_graphics/big_barrel/circle") {

  auto scanner = BarrelBuilder::make_big_barrel();

  std::ofstream out("test_output/graph_circle.m");
  if (!out) {
    FAIL("cannot open file");
  }

  MathematicaGraphics graphics(out);

  graphics.add(scanner);
  graphics.add_circle(0.400);

  for (Detector& d : scanner) {
    auto center = d.center();
    graphics.add_circle(center, 0.015);
  }
}

TEST("common/mathematica_graphics/big_barrel/pixel") {

  auto scanner = BarrelBuilder::make_big_barrel();
  using Point = PET2D::Point<F>;

  std::ofstream out("test_output/graph_pixels.m");
  if (!out) {
    FAIL("cannot open file");
  }

  MathematicaGraphics graphics(out);

  graphics.add(scanner);
  const int n_columns = 20;
  const int n_rows = 20;
  PET2D::PixelGrid<F, S> grid(n_columns, n_rows, 0.01, Point(-0.1, -0.1));

  for (int ix = 0; ix < n_columns; ++ix) {
    for (int iy = 0; iy < n_rows; ++iy) {
      graphics.add_pixel(grid, ix, iy);
    }
  }
}
