#include "util/test.h"

#include "scanner_builder.h"
#include "generic_scanner.h"
#include "square_detector.h"

using Builder = PET2D::Barrel::ScannerBuilder<
    PET2D::Barrel::
        GenericScanner<PET2D::Barrel::SquareDetector<float>, short, 128>>;

TEST("2d/barrel/scanner_builder/single_ring/symmetry") {

  auto detector = Builder::build_single_ring(200.0, 8, 0.007, 0.019);

  auto symmetry_descriptor = detector.symmetry_descriptor();

  // first ring
  for (short d = 0; d < 8; d++)
    REQUIRE(symmetry_descriptor->symmetric_detector(d, 0) == d);

  REQUIRE(symmetry_descriptor->symmetric_detector(1, 1) == 7);
  REQUIRE(symmetry_descriptor->symmetric_detector(1, 2) == 3);
  REQUIRE(symmetry_descriptor->symmetric_detector(1, 3) == 5);
  REQUIRE(symmetry_descriptor->symmetric_detector(1, 4) == 1);
  REQUIRE(symmetry_descriptor->symmetric_detector(1, 5) == 3);
  REQUIRE(symmetry_descriptor->symmetric_detector(1, 6) == 7);
  REQUIRE(symmetry_descriptor->symmetric_detector(1, 7) == 5);

  REQUIRE(symmetry_descriptor->symmetric_detector(6, 1) == 2);
  REQUIRE(symmetry_descriptor->symmetric_detector(6, 2) == 6);
  REQUIRE(symmetry_descriptor->symmetric_detector(6, 3) == 2);
  REQUIRE(symmetry_descriptor->symmetric_detector(6, 4) == 4);
  REQUIRE(symmetry_descriptor->symmetric_detector(6, 5) == 0);
  REQUIRE(symmetry_descriptor->symmetric_detector(6, 6) == 4);
  REQUIRE(symmetry_descriptor->symmetric_detector(6, 7) == 0);
}

TEST("2d/barrel/scanner_builder/multi_ring/symmetry") {

  auto detector = Builder::build_multiple_rings(
      { 425, 475, 525 }, { 0.0, 0.5, 0.0 }, { 8, 12, 24 }, 0.007, 0.019);

  auto symmetry_descriptor = detector.symmetry_descriptor();

  // first ring
  for (short d = 0; d < 8; d++)
    REQUIRE(symmetry_descriptor->symmetric_detector(d, 0) == d);

  REQUIRE(symmetry_descriptor->symmetric_detector(1, 1) == 7);
  REQUIRE(symmetry_descriptor->symmetric_detector(1, 2) == 3);
  REQUIRE(symmetry_descriptor->symmetric_detector(1, 3) == 5);
  REQUIRE(symmetry_descriptor->symmetric_detector(1, 4) == 1);
  REQUIRE(symmetry_descriptor->symmetric_detector(1, 5) == 3);
  REQUIRE(symmetry_descriptor->symmetric_detector(1, 6) == 7);
  REQUIRE(symmetry_descriptor->symmetric_detector(1, 7) == 5);

  REQUIRE(symmetry_descriptor->symmetric_detector(6, 1) == 2);
  REQUIRE(symmetry_descriptor->symmetric_detector(6, 2) == 6);
  REQUIRE(symmetry_descriptor->symmetric_detector(6, 3) == 2);
  REQUIRE(symmetry_descriptor->symmetric_detector(6, 4) == 4);
  REQUIRE(symmetry_descriptor->symmetric_detector(6, 5) == 0);
  REQUIRE(symmetry_descriptor->symmetric_detector(6, 6) == 4);
  REQUIRE(symmetry_descriptor->symmetric_detector(6, 7) == 0);

  // second ring
  for (short d = 8; d < 8 + 12; d++)
    REQUIRE(symmetry_descriptor->symmetric_detector(d, 0) == d);

  REQUIRE(symmetry_descriptor->symmetric_detector(8, 1) == 19);
  REQUIRE(symmetry_descriptor->symmetric_detector(8, 2) == 13);
  REQUIRE(symmetry_descriptor->symmetric_detector(8, 3) == 14);
  REQUIRE(symmetry_descriptor->symmetric_detector(8, 4) == 10);
  REQUIRE(symmetry_descriptor->symmetric_detector(8, 5) == 11);
  REQUIRE(symmetry_descriptor->symmetric_detector(8, 6) == 17);
  REQUIRE(symmetry_descriptor->symmetric_detector(8, 7) == 16);

  REQUIRE(symmetry_descriptor->symmetric_detector(15, 1) == 12);
  REQUIRE(symmetry_descriptor->symmetric_detector(15, 2) == 18);
  REQUIRE(symmetry_descriptor->symmetric_detector(15, 3) == 9);
  REQUIRE(symmetry_descriptor->symmetric_detector(15, 4) == 15);
  REQUIRE(symmetry_descriptor->symmetric_detector(15, 5) == 18);
  REQUIRE(symmetry_descriptor->symmetric_detector(15, 6) == 12);
  REQUIRE(symmetry_descriptor->symmetric_detector(15, 7) == 9);

  // third ring
  for (short d = 20; d < 8 + 12 + 24; d++)
    REQUIRE(symmetry_descriptor->symmetric_detector(d, 0) == d);

  REQUIRE(symmetry_descriptor->symmetric_detector(20, 1) == 20);
  REQUIRE(symmetry_descriptor->symmetric_detector(20, 2) == 32);
  REQUIRE(symmetry_descriptor->symmetric_detector(20, 3) == 32);
  REQUIRE(symmetry_descriptor->symmetric_detector(20, 4) == 26);
  REQUIRE(symmetry_descriptor->symmetric_detector(20, 5) == 26);
  REQUIRE(symmetry_descriptor->symmetric_detector(20, 6) == 38);
  REQUIRE(symmetry_descriptor->symmetric_detector(20, 7) == 38);
}
