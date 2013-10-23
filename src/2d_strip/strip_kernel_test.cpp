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

TEST_CASE("sensitivity tests/square", "square") {

  StripDetector<double> detector(500, 1000, 200, 200, 5, 5, 10, 63);

  CHECK(detector.sensitivity(0.0, 0.0) == Approx(0.5).epsilon(1e-13));
  CHECK(detector.sensitivity(0.0, 50.0) ==
        Approx(0.46652458328685176).epsilon(1e-13));
  CHECK(detector.sensitivity(100.0, -50.0) ==
        Approx(0.4410019151324715).epsilon(1e-13));
  CHECK(detector.sensitivity(-200, -450) ==
        Approx(0.07526632771111386).epsilon(1e-13));
}

TEST_CASE("sensitivity tests/non square", "non square") {

  StripDetector<double> detector(450, 200, 200, 200, 5, 5, 10, 63);

  CHECK(detector.sensitivity(0.0, 0.0) ==
        Approx(0.1392089745461279).epsilon(1e-13));
  CHECK(detector.sensitivity(0.0, 50.0) ==
        Approx(0.07044657495455454).epsilon(1e-13));
  CHECK(detector.sensitivity(100.0, -50.0) ==
        Approx(0.07402517367717103).epsilon(1e-13));
  CHECK(detector.sensitivity(-200, -70) ==
        Approx(0.05269621503719814).epsilon(1e-13));
}

#if DONT_USE
TEST_CASE("kernel tests constructor 1", "kernel") {
  StripDetector<double> detector(500, 1000, 200, 200, 5, 5, 10, 63);

  Reconstruction<double> reconstructor(1, detector);

  check(1.1372205719261035e-7, 0.0, 0.0, 0.0, 0.0, reconstructor);

  check(1.99620227633633e-8, 0.0, 0.0, 10.0, 13.0, reconstructor);

  check(5.5729829923449995e-8, 100.0, 45.0 * degree, 0.0, 0.0, reconstructor);

  check(3.12537857516921e-11, 100.0, 45.0 * degree, -20.0, 7.0, reconstructor);

  check(7.993589560016591e-8, -10.0, -13.0 * degree, -2.0, -5.0, reconstructor);
}
#endif

TEST_CASE("kernel tests constructor 2", "kernel") {

  Reconstruction<double> reconstructor(1, 500, 100, 200, 5, 10, 63);

  check(1.1372205719261035e-7, 0.0, 0.0, 0.0, 0.0, reconstructor);

  check(1.99620227633633e-8, 0.0, 0.0, 10.0, 13.0, reconstructor);

  check(5.5729829923449995e-8, 100.0, 45.0 * degree, 0.0, 0.0, reconstructor);

  check(3.12537857516921e-11, 100.0, 45.0 * degree, -20.0, 7.0, reconstructor);

  check(7.993589560016591e-8, -10.0, -13.0 * degree, -2.0, -5.0, reconstructor);
}
