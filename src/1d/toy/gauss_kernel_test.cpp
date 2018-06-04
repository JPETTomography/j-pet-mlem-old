#include "util/test.h"

#include "1d/toy/gauss_kernel.h"

using Kernel = PET1D::Toy::GaussKernel<float>;

TEST("1d/toy/gauss_kernel") {

  Kernel kernel(0.75);

  REQUIRE(kernel(0, 1) == Approx(0.218680099568));
  REQUIRE(kernel(0.53, 2) == Approx(0.07792125911126863));
}
