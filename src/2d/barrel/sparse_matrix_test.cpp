#include <iostream>

#include "util/test.h"

#include "sparse_matrix.h"
#include "2d/geometry/pixel.h"
#include "2d/barrel/lor.h"

TEST("SparseMatrix/symmetric_lor") {

  PET2D::Barrel::SparseMatrix<PET2D::Pixel<short>,
                              PET2D::Barrel::LOR<short>,
                              short,
                              int> matrix(128, 24, 0, 1);

  short detector = 1;
  REQUIRE(matrix.symmetric_detector(detector, 0) == 1);
  REQUIRE(matrix.symmetric_detector(detector, 1) == 23);
  REQUIRE(matrix.symmetric_detector(detector, 2) == 11);
  REQUIRE(matrix.symmetric_detector(detector, 3) == 13);
  REQUIRE(matrix.symmetric_detector(detector, 4) == 5);
  REQUIRE(matrix.symmetric_detector(detector, 5) == 7);
  REQUIRE(matrix.symmetric_detector(detector, 6) == 19);
  REQUIRE(matrix.symmetric_detector(detector, 7) == 17);
}
