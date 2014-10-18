#include "catch.hpp"

#include "matrix.h"

using namespace Math;

TEST_CASE("math/matrix/ctor") {
  Matrix<3> m{ 1, 2, 3,  //
               3, 4, 5,  //
               5, 6, 7 };
}
