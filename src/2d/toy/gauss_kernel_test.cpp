#include "util/test.h"

#include "2d/toy/gauss_kernel.h"

TEST("2d/toy/gauss_kernel/value") {

  PET2D::Toy::GaussKernel<float> kernel(1.2, 2.7);

  REQUIRE(1.2 == Approx(kernel.sigma_x));
  REQUIRE(2.7 == Approx(kernel.sigma_y));

  auto res = kernel(0.5, 0.7);

  REQUIRE(res == Approx(0.043549214).epsilon(1e-7));
}

TEST("2d/toy/gauss_kernel/bb") {

  PET2D::Toy::GaussKernel<float> kernel(1.2, 2.7);

  REQUIRE(1.2 == Approx(kernel.sigma_x));
  REQUIRE(2.7 == Approx(kernel.sigma_y));

  REQUIRE(kernel.three_sigma_bb_x() == Approx(3 * 1.2).epsilon(1e-7));
  REQUIRE(kernel.three_sigma_bb_y() == Approx(3 * 2.7).epsilon(1e-7));
}
