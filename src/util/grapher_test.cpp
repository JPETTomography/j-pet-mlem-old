#include <fstream>

#include "util/test.h"

#include "grapher.h"

#include "2d/barrel/square_detector.h"
#include "2d/barrel/barrel_builder.h"
#include "2d/geometry/line_segment.h"

TEST("grapher/detector") {
  using Detector = PET2D::Barrel::SquareDetector<float>;
  using F = Detector::F;

  std::ofstream out("test_output/graph_detector.m");
  if (!out) {
    FAIL("cannot open file");
  }

  Graphics<F> graphics(out);

  Detector detector(0.007, 0.019, 0);

  graphics.add(detector);
}

TEST("grapher/BigBarel") {

  PET2D::Barrel::BigBarrelType scanner = PET2D::Barrel::buildBigBarrel();
  using F = PET2D::Barrel::BigBarrelType::F;

  std::ofstream out("test_output/graph_scanner.m");
  if (!out) {
    FAIL("cannot open file");
  }

  Graphics<F> graphics(out);

  graphics.add(scanner);
}

TEST("grapher/BigBarel/lor") {

  PET2D::Barrel::BigBarrelType scanner = PET2D::Barrel::buildBigBarrel();
  using F = PET2D::Barrel::BigBarrelType::F;

  std::ofstream out("test_output/graph_lor.m");
  if (!out) {
    FAIL("cannot open file");
  }

  Graphics<F> graphics(out);

  graphics.add(scanner);
  graphics.add(scanner, PET2D::Barrel::LOR<int>(65, 0));
}

TEST("grapher/BigBarel/segment") {

  PET2D::Barrel::BigBarrelType scanner = PET2D::Barrel::buildBigBarrel();
  using F = PET2D::Barrel::BigBarrelType::F;
  using Point = PET2D::Point<F>;

  std::ofstream out("test_output/graph_segment.m");
  if (!out) {
    FAIL("cannot open file");
  }

  Graphics<F> graphics(out);

  graphics.add(scanner);

  PET2D::LineSegment<F> segment(Point(-0.400, 0), Point(0, 0.400));

  graphics.add(segment);
}

TEST("grapher/BigBarel/circle") {

  PET2D::Barrel::BigBarrelType scanner = PET2D::Barrel::buildBigBarrel();
  using Detector = PET2D::Barrel::BigBarrelType::Detector;
  using F = PET2D::Barrel::BigBarrelType::F;

  std::ofstream out("test_output/graph_circle.m");
  if (!out) {
    FAIL("cannot open file");
  }

  Graphics<F> graphics(out);

  graphics.add(scanner);
  graphics.addCircle(0.400);

  for (Detector& d : scanner) {
    auto center = d.center();
    graphics.addCircle(center, 0.015);
  }
}
