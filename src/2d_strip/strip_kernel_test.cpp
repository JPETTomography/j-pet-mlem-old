#include "catch.hpp"
#include <vector>
#include <cmath>

#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"

#include "reconstruction.h"
#include "kernel.h"

const double degree = M_PI / 180.0;

template <typename F>
void check(double ref,
           double y,
           double angle,  // radians
           double dy,
           double dz,
           StripDetector<F>& detector) {

  Kernel<F> kernel = Kernel<F>();

  double tan_value = tan(angle);
  double inv_cos = 1.0 / cos(angle);
  double inv_cos_sq = inv_cos * inv_cos;
  Point<> delta(dy, dz);

  double value = kernel(y,
                        tan_value,
                        inv_cos,
                        inv_cos_sq,
                        detector.radius,
                        delta,
                        detector.inv_cor_mat_diag,
                        detector.sqrt_det_cor_mat());
  CHECK(value == Approx(ref).epsilon(1e-13));
}

TEST_CASE("strip/sensitivity/square") {

  StripDetector<> detector(500, 1000, 200, 200, 5, 5, 10, 63);

  CHECK(detector.sensitivity(Point<>(0.0, 0.0)) == Approx(0.5).epsilon(1e-13));
  CHECK(detector.sensitivity(Point<>(0.0, 50.0)) ==
        Approx(0.46652458328685176).epsilon(1e-13));
  CHECK(detector.sensitivity(Point<>(100.0, -50.0)) ==
        Approx(0.4410019151324715).epsilon(1e-13));
  CHECK(detector.sensitivity(Point<>(-200, -450)) ==
        Approx(0.07526632771111386).epsilon(1e-13));
}

TEST_CASE("strip/sensitivity/non_square") {

  StripDetector<> detector(450, 200, 200, 200, 5, 5, 10, 63);

  CHECK(detector.sensitivity(Point<>(0.0, 0.0)) ==
        Approx(0.1392089745461279).epsilon(1e-13));
  CHECK(detector.sensitivity(Point<>(0.0, 50.0)) ==
        Approx(0.07044657495455454).epsilon(1e-13));
  CHECK(detector.sensitivity(Point<>(100.0, -50.0)) ==
        Approx(0.07402517367717103).epsilon(1e-13));
  CHECK(detector.sensitivity(Point<>(-200, -70)) ==
        Approx(0.05269621503719814).epsilon(1e-13));
}

#if DONT_TEST
TEST_CASE("strip/kernel/ctor1", "[ctor]") {

  StripDetector<> detector(500, 1000, 200, 200, 5, 5, 10, 63);

  check(1.1372205719261035e-7, 0.0, 0.0, 0.0, 0.0, detector);
  check(1.99620227633633e-8, 0.0, 0.0, 10.0, 13.0, detector);
  check(5.5729829923449995e-8, 100.0, 45.0 * degree, 0.0, 0.0, detector);
  check(3.12537857516921e-11, 100.0, 45.0 * degree, -20.0, 7.0, detector);
  check(7.993589560016591e-8, -10.0, -13.0 * degree, -2.0, -5.0, detector);
}
#endif

TEST_CASE("strip/kernel/ctor2", "[ctor]") {

  StripDetector<> detector(500, 1000, 200, 200, 5, 5, 10, 63);

#if DONT_TEST
  check(1.1372205719261035e-7, 0.0, 0.0, 0.0, 0.0, detector);
#endif
  check(1.99620227633633e-8, 0.0, 0.0, 10.0, 13.0, detector);
#if DONT_TEST
  check(5.5729829923449995e-8, 100.0, 45.0 * degree, 0.0, 0.0, detector);
  check(3.12537857516921e-11, 100.0, 45.0 * degree, -20.0, 7.0, detector);
  check(7.993589560016591e-8, -10.0, -13.0 * degree, -2.0, -5.0, detector);
#endif
}
