#include "catch.hpp"

#include "../point.h"

#include"tor.h"

typedef double FLOAT;

TEST_CASE("tor/create", "TOR create") {
  Point<FLOAT> c1(450, 0);
  Point<FLOAT> c2(-450, 0);

  ToR<FLOAT> tor(c1, 0.0, 5, 19, c2, M_PI, 5, 19);

}
