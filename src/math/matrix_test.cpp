#include "catch.hpp"

#include "matrix.h"

TEST_CASE("matrix/init", "matrix initialization") {
  Matrix<3> m{ 1, 2, 3,  //
               3, 4, 5,  //
               5, 6, 7 };
}
