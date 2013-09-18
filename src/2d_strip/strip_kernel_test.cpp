#include "catch.hpp"
#include <vector>
#include <cmath>

#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"

#include"reconstruction.h"

typedef Reconstruction<double>::Point Point;

const double degree = M_PI / 180.0;

void check(double ref,
           double y,
           double angle,  // radians
           double dy,
           double dz,
           Reconstruction<double>& rec) {

  double tan_value = tan(angle);
  double inv_cos = 1.0 / cos(angle);
  double inv_cos_sq = inv_cos * inv_cos;
  Point delta(dy, dz);
  double value = rec.kernel(y, tan_value, inv_cos, inv_cos_sq, delta);
  CHECK(value == Approx(ref).epsilon(1e-13));
}

TEST_CASE("kernel tests", "kernel") {

  Reconstruction<double> reconstructor(1, 500, 1000, 200, 5, 10, 63);

  check(1.1372205719261035e-7, 0.0, 0.0, 0.0, 0.0, reconstructor);

  check(1.99620227633633e-8, 0.0, 0.0, 10.0, 13.0, reconstructor);

  check(5.5729829923449995e-8, 100.0, 45.0 * degree, 0.0, 0.0, reconstructor);

  check(3.12537857516921e-11, 100.0, 45.0 * degree, -20.0, 7.0, reconstructor);

  check(7.993589560016591e-8, -10.0, -13.0 * degree, -2.0, -5.0, reconstructor);
}
