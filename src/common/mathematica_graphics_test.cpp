#include <sstream>

#include "util/test.h"

#include "mathematica_graphics.h"

#include "2d/barrel/square_detector.h"
#include "2d/barrel/scanner_builder.h"
#include "2d/geometry/line_segment.h"
#include "2d/geometry/pixel_grid.h"
#include "2d/barrel/options.h"

#include "common/types.h"

using Detector = PET2D::Barrel::SquareDetector<F>;
using Scanner = PET2D::Barrel::GenericScanner<Detector, S>;
using ScannerBuilder = PET2D::Barrel::ScannerBuilder<Scanner>;
using MathematicaGraphics = Common::MathematicaGraphics<F>;

static Scanner make_big_barrel() {
  std::ifstream in_cfg("samples/config/big.cfg");
  cmdline::parser cl;
  PET2D::Barrel::add_scanner_options(cl);
  cl.parse(in_cfg);
  return ScannerBuilder::build_multiple_rings(
      PET2D_BARREL_SCANNER_CL(cl, Detector::F));
}

TEST("common/mathematica_graphics/detector") {
  std::stringstream out;
  Detector detector(0.007, 0.019, 0);
  {
    MathematicaGraphics graphics(out);
    graphics.add(detector);
  }
  REQUIRE(out.str() ==
          "Graphics[{\n"
          "{Polygon[{\n"
          "  {0.00350000010803, 0.00949999969453},\n"
          "  {0.00350000010803, -0.00949999969453},\n"
          "  {-0.00350000010803, -0.00949999969453},\n"
          "  {-0.00350000010803, 0.00949999969453}}]}}]\n");
}

#define TEST_LINE(out, line, text) \
  std::getline(out, line);         \
  REQUIRE(line == text);

TEST("common/mathematica_graphics/big_barrel") {
  std::stringstream out;
  auto scanner = make_big_barrel();
  {
    MathematicaGraphics graphics(out);
    graphics.add(scanner);
  }
  std::string line;
  TEST_LINE(out, line, "Graphics[{");
  TEST_LINE(out, line, "{{Polygon[{");
  TEST_LINE(out, line, "  {0.439500004053, -0.00350001943298},");
}

TEST("common/mathematica_graphics/big_barrel/lor") {
  std::stringstream out;
  auto scanner = make_big_barrel();
  {
    MathematicaGraphics graphics(out);
    graphics.add(scanner);
    graphics.add(scanner, PET2D::Barrel::LOR<int>(65, 0));
  }
  std::string line;
  TEST_LINE(out, line, "Graphics[{");
  TEST_LINE(out, line, "{{Polygon[{");
  TEST_LINE(out, line, "  {0.439500004053, -0.00350001943298},");
}

TEST("common/mathematica_graphics/big_barrel/segment") {
  std::stringstream out;
  auto scanner = make_big_barrel();
  using Point = PET2D::Point<F>;
  {
    MathematicaGraphics graphics(out);
    graphics.add(scanner);
    PET2D::LineSegment<F> segment(Point(-0.400, 0), Point(0, 0.400));
    graphics.add(segment);
  }
  std::string line;
  TEST_LINE(out, line, "Graphics[{");
  TEST_LINE(out, line, "{{Polygon[{");
  TEST_LINE(out, line, "  {0.439500004053, -0.00350001943298},");
}

TEST("common/mathematica_graphics/big_barrel/circle") {
  std::stringstream out;
  auto scanner = make_big_barrel();
  {
    MathematicaGraphics graphics(out);
    graphics.add(scanner);
    graphics.add_circle(0.400);
    for (Detector& d : scanner) {
      auto center = d.center();
      graphics.add_circle(center, 0.015);
    }
  }
  std::string line;
  TEST_LINE(out, line, "Graphics[{");
  TEST_LINE(out, line, "{{Polygon[{");
  TEST_LINE(out, line, "  {0.439500004053, -0.00350001943298},");
}

TEST("common/mathematica_graphics/big_barrel/pixel") {
  std::stringstream out;
  auto scanner = make_big_barrel();
  using Point = PET2D::Point<F>;
  {
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
  std::string line;
  TEST_LINE(out, line, "Graphics[{");
  TEST_LINE(out, line, "{{Polygon[{");
  TEST_LINE(out, line, "  {0.439500004053, -0.00350001943298},");
}
