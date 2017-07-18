#include <vector>

#include "util/test.h"
#include "util/array.h"

#include "symmetry_descriptor.h"
#include "2d/geometry/point.h"
#include "common/types.h"

TEST("Symmetry transformations") {
  using Point = PET2D::Point<F>;
  using Transformation = PET2D::Transformation<F>;
  using SymmetryDescriptor = PET2D::Barrel::SymmetryDescriptor<S>;

  std::vector<Transformation> tr;
  tr.reserve(SymmetryDescriptor::EIGHT);
  for (int i = 0; i < SymmetryDescriptor::EIGHT; i++) {
    new (&tr[i])
        Transformation(SymmetryDescriptor::symmetry_transformation<F>(i));
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
