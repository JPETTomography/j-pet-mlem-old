#include "catch.hpp"
#include <vector>
#include <cmath>

#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"

#include"reconstruction.h"

typedef Reconstruction<float>::Point Point;

const float degree = M_PI / 180.0;
void check(float ref,
           float y,
           float angle,  // radians
           float dy,
           float dz,
           Reconstruction<float>& rec) {

  float tan_value = tan(angle);
  float inv_cos = 1.0 / cos(angle);
  float inv_cos_sq = inv_cos * inv_cos;
  Point delta(dy, dz);
  float value = rec.kernel(y, tan_value, inv_cos, inv_cos_sq, delta);
  CHECK(value == Approx(ref).epsilon(1e-8));
}

TEST_CASE("kernel tests", "kernel") {

  Reconstruction<float> reconstructor(1, 500, 1000, 200, 5, 10, 63);

  check(1.1372205719261035e-7, 0.0, 0.0, 0.0, 0.0, reconstructor);

  check(1.99620227633633e-8, 0.0, 0.0, 10.0, 13.0, reconstructor);

  check(5.5729829923449995e-8, 100.0, 45.0 * degree, 0.0, 0.0, reconstructor);

  check(3.12537857516921e-11, 100.0, 45.0 * degree, -20.0, 7.0, reconstructor);

  check(7.993589560016591e-8, -10.0, -13.0 * degree, -2.0, -5.0, reconstructor);
}
