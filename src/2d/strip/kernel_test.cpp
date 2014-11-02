#include <vector>
#include <cmath>

#include "util/test.h"
#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"

#include "reconstruction.h"
#include "kernel.h"

using namespace PET2D;
using namespace PET2D::Strip;

const double degree = M_PI / 180.0;

template <typename F>
void check(F ref, F y, F angle, F dy, F dz, F R, const Kernel<F>& kernel) {
  F tangent = std::tan(angle);
  F secant = 1 / std::cos(angle);
  Point<F> delta(dy, dz);
  CHECK(kernel(y, tangent, secant, R, delta) == Approx(ref).epsilon(1e-13));
}

TEST_CASE("strip/sensitivity/square") {

  Detector<> detector(500, 1000, 200, 200, 5, 5, 10, 63);

  CHECK(detector.sensitivity(Point<>(0.0, 0.0)) == 0.5_e13);
  CHECK(detector.sensitivity(Point<>(0.0, 50.0)) == 0.46652458328685176_e13);
  CHECK(detector.sensitivity(Point<>(100.0, -50.0)) == 0.4410019151324715_e13);
  CHECK(detector.sensitivity(Point<>(-200, -450)) == 0.07526632771111386_e13);
}

TEST_CASE("strip/sensitivity/non_square") {

  Detector<> detector(450, 200, 200, 200, 5, 5, 10, 63);

  CHECK(detector.sensitivity(Point<>(0.0, 0.0)) == 0.1392089745461279_e13);
  CHECK(detector.sensitivity(Point<>(0.0, 50.0)) == 0.07044657495455454_e13);
  CHECK(detector.sensitivity(Point<>(100.0, -50.0)) == 0.07402517367717103_e13);
  CHECK(detector.sensitivity(Point<>(-200, -70)) == 0.05269621503719814_e13);
}

#if DONT_TEST
TEST_CASE("strip/kernel/ctor1") {

  Detector<> detector(500, 1000, 200, 200, 5, 5, 10, 63);
  Kernel<> kernel(detector.sigma_z, detector.sigma_dl);
  double R = detector.radius;

  check(1.1372205719261035e-7, 0.0, 0.0, 0.0, 0.0, R, kernel);
  check(1.99620227633633e-8, 0.0, 0.0, 10.0, 13.0, R, kernel);
  check(5.5729829923449995e-8, 100.0, 45.0 * degree, 0.0, 0.0, R, kernel);
  check(3.12537857516921e-11, 100.0, 45.0 * degree, -20.0, 7.0, R, kernel);
  check(7.993589560016591e-8, -10.0, -13.0 * degree, -2.0, -5.0, R, kernel);
}
#endif

TEST_CASE("strip/kernel/ctor2") {

  Detector<> detector(500, 1000, 200, 200, 5, 5, 10, 63);
  Kernel<> kernel(detector.sigma_z, detector.sigma_dl);
  double R = detector.radius;

#if DONT_TEST
  check(1.1372205719261035e-7, 0.0, 0.0, 0.0, 0.0, R, kernel);
#endif
  check(1.99620227633633e-8, 0.0, 0.0, 10.0, 13.0, R, kernel);
#if DONT_TEST
  check(5.5729829923449995e-8, 100.0, 45.0 * degree, 0.0, 0.0, R, kernel);
  check(3.12537857516921e-11, 100.0, 45.0 * degree, -20.0, 7.0, R, kernel);
  check(7.993589560016591e-8, -10.0, -13.0 * degree, -2.0, -5.0, R, kernel);
#endif
}

TEST_CASE("strip/kernel/bbox") {

  Detector<> detector(500, 1000, 200, 200, 5, 5, 10, 63);
  Kernel<> kernel(detector.sigma_z, detector.sigma_dl);
  double R = detector.radius;

  struct {
    double angle;
    double bby_value;
    double bbz_value;
  } v[] = { { 0.0470448, 94.3954, 21.6737 },
            { -0.594145, 78.3053, 56.9959 },
            { 0.20029, 92.6108, 28.3458 },
            { -0.571667, 79.4745, 55.3539 },
            { -0.420542, 86.266, 44.0276 } };

  for (size_t i = 0; i < sizeof(v) / sizeof(*v); ++i) {
    auto inv_pow_sigma_dl = (1.0 / (detector.sigma_dl * detector.sigma_dl));
    auto inv_pow_sigma_z = (1.0 / (detector.sigma_z * detector.sigma_z));
    auto angle = v[i].angle;
    auto bby_value = v[i].bby_value;
    auto bbz_value = v[i].bbz_value;

    auto A = ((4.0 / (cos(angle) * cos(angle))) * inv_pow_sigma_dl) +
             (2.0 * tan(angle) * tan(angle) * inv_pow_sigma_z);
    auto B = -4.0 * tan(angle) * inv_pow_sigma_z;
    auto C = 2.0 * inv_pow_sigma_z;
    auto B_2 = (B / 2.0) * (B / 2.0);

    CHECK(kernel.bb_y(A, C, B_2) == Approx(bby_value).epsilon(1e-4));
    CHECK(kernel.bb_z(A, C, B_2) == Approx(bbz_value).epsilon(1e-4));
  }
}
