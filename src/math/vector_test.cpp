#include "catch.hpp"

#include "vector.h"

TEST_CASE("vector/ctor", "[ctor]") {
  Vector<3> v{ 1, 2, 3 };
}
