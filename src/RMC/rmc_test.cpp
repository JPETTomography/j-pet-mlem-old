#include "catch.hpp"

#include "../point.h"

#include "tor.h"

typedef double FLOAT;

TEST_CASE("tor/create", "TOR create") {
  Point<FLOAT> c[2] = { { 450, 0 }, { -450, 0 } };
  FLOAT a[2] = { 0.0, M_PI };
  FLOAT w[2] = { 5, 6 };
  FLOAT h[2] = { 19, 20 };

  ToR<FLOAT> tor(c[0], 0.0, 5, 19, c[1], M_PI, 6, 20);

  for (int i = 0; i < 2; ++i) {
    REQUIRE(tor.center(i).x == Approx(c[i].x));
    REQUIRE(tor.center(i).y == Approx(c[i].y));
    REQUIRE(tor.angle(i) == Approx(a[i]));
    REQUIRE(tor.width(i) == Approx(w[i]));
    REQUIRE(tor.height(i) == Approx(h[i]));
  }

}
