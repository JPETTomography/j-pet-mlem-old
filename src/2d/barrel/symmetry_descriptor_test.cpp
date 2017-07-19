#include <vector>

#include "util/test.h"
#include "util/array.h"

#include "symmetry_descriptor.h"
#include "2d/geometry/point.h"

#include "scanner_builder.h"
#include "generic_scanner.h"
#include "square_detector.h"
#include "2d/geometry/find_symmetry.h"

#include "common/types.h"

using Point = PET2D::Point<F>;
using Transformation = PET2D::Transformation<F>;
using SymmetryDescriptor = PET2D::Barrel::SymmetryDescriptor<S>;

TEST("Symmetry transformations") {

  std::vector<Transformation> tr;
  tr.reserve(SymmetryDescriptor::EIGHT);
  for (short i = 0; i < SymmetryDescriptor::EIGHT; i++) {
    new (&tr[i]) Transformation(symmetry_transformation<F>(i));
  }

  Point p(0.3, 0.7);

  REQUIRE(tr[0](p) == p);
  REQUIRE(tr[1](p) == Point(-p.x, p.y));
  REQUIRE(tr[2](p).approx_equal(Point(p.x, -p.y)));
  REQUIRE(tr[3](p).approx_equal(Point(-p.x, -p.y)));

  REQUIRE(tr[4](p).approx_equal(Point(p.y, p.x)));
  REQUIRE(tr[5](p).approx_equal(Point(p.y, -p.x)));
  REQUIRE(tr[6](p).approx_equal(Point(-p.y, p.x)));
  REQUIRE(tr[7](p).approx_equal(Point(-p.y, -p.x)));
}

TEST("Find symmetry") {
  using Builder = PET2D::Barrel::ScannerBuilder<
      PET2D::Barrel::GenericScanner<PET2D::Barrel::SquareDetector<F>, S, 128>>;

  auto detector = Builder::build_single_ring(200.0, 8, 0.007, 0.019);

  auto symmetry_descriptor = detector.symmetry_descriptor();

  // first ring
  for (short d = 0; d < 8; d++)
    REQUIRE(symmetry_descriptor.symmetric_detector(d, 0) == d);

  REQUIRE(symmetry_descriptor.symmetric_detector(1, 1) == 7);
  REQUIRE(symmetry_descriptor.symmetric_detector(1, 2) == 3);
  REQUIRE(symmetry_descriptor.symmetric_detector(1, 3) == 5);
  REQUIRE(symmetry_descriptor.symmetric_detector(1, 4) == 1);
  REQUIRE(symmetry_descriptor.symmetric_detector(1, 5) == 3);
  REQUIRE(symmetry_descriptor.symmetric_detector(1, 6) == 7);
  REQUIRE(symmetry_descriptor.symmetric_detector(1, 7) == 5);

  REQUIRE(symmetry_descriptor.symmetric_detector(6, 1) == 2);
  REQUIRE(symmetry_descriptor.symmetric_detector(6, 2) == 6);
  REQUIRE(symmetry_descriptor.symmetric_detector(6, 3) == 2);
  REQUIRE(symmetry_descriptor.symmetric_detector(6, 4) == 4);
  REQUIRE(symmetry_descriptor.symmetric_detector(6, 5) == 0);
  REQUIRE(symmetry_descriptor.symmetric_detector(6, 6) == 4);
  REQUIRE(symmetry_descriptor.symmetric_detector(6, 7) == 0);

  for (S s = 0; s < SymmetryDescriptor::EIGHT; s++) {
    for (S d = 0; d < detector.size(); d++) {
      REQUIRE(find_symmetric(detector, s, d) ==
              symmetry_descriptor.symmetric_detector(s, d));
    }
  }
}
