#include <fstream>

#include "util/test.h"

#include "grapher.h"

#include "2d/barrel/square_detector.h"
#include "2d/barrel/barrel_builder.h"

TEST("grapher/detector") {
  using Detector = PET2D::Barrel::SquareDetector<float>;
  using F = Detector::F;

  std::ofstream out("test_output/graph.m");
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

  std::ofstream out("test_output/graph.m");
  if (!out) {
    FAIL("cannot open file");
  }

  Graphics<F> graphics(out);

  graphics.add(scanner);
}
