#include <iostream>

#include "util/test.h"

#include "sparse_matrix.h"
#include "2d/geometry/pixel.h"
#include "2d/barrel/lor.h"

TEST("2d/barrel/sparse_matrix/symmetric_lor") {

  PET2D::Barrel::SparseMatrix<PET2D::Pixel<short>,
                              PET2D::Barrel::LOR<short>,
                              short,
                              int> matrix(128, 24, 0, 1);

  short detector = 1;
}
