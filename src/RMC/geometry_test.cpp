#include "catch.hpp"

#include "../geometry/geometry.h"

using namespace Geometry;

TEST_CASE("geometry/vector", "Vector") {
  double init[2] = { 0.0, 1.0 };
  Vector<2> v(init);

  for (int i = 0; i < 2; ++i)
    REQUIRE(v[i] == Approx(init[i]));
}
