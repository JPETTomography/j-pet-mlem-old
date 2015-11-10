#include <cmath>

#include "util/test.h"

#include "common/types.h"
#include "2d/strip/gausian_kernel.h"

using Kernel = PET2D::Strip::GaussianKernel<F>;
using Vector = PET2D::Vector<F>;

const float sigma_z = 0.01;
const float sigma_dl = 0.04;
const float epsilon = 1e-6;

TEST("2d/strip/gaussian_kernel") {

  Kernel gaussian(sigma_z, sigma_dl);
  F R = 0.35;

  SECTION("event 1") {
    F angle = M_PI / 6;
    F sec = 1 / std::cos(angle);
    F tan = std::tan(angle);

    F y = 0.1;

    CHECK(gaussian(y, tan, sec, R, Vector(0, 0)) ==
          Approx(4.6076243632716500e07).epsilon(epsilon));
    CHECK(gaussian(y, tan, sec, R, Vector(0.05, 0.1)) ==
          Approx(1.3180240903945950).epsilon(epsilon));
#ifndef __clang__
    // FIXME: Clang provides invalid value of 8461.84668f here
    CHECK(gaussian(y, tan, sec, R, Vector(-0.05, -0.05)) ==
          Approx(8461.8357488172500000).epsilon(epsilon));
#endif
  }

  SECTION("event 2") {
    F angle = M_PI / 4;
    F sec = 1 / std::cos(angle);
    F tan = std::tan(angle);
    F y = -0.3;

    CHECK(gaussian(y, tan, sec, R, Vector(0, 0)) ==
          Approx(2.3969145992977610e07).epsilon(epsilon));
    CHECK(gaussian(y, tan, sec, R, Vector(0.05, 0.1)) ==
          Approx(0.0).epsilon(epsilon));
    CHECK(gaussian(y, tan, sec, R, Vector(-0.05, -0.05)) ==
          Approx(42993.19771499309000).epsilon(epsilon));
  }
}
