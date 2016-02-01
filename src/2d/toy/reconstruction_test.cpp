#include <sstream>

#include "util/test.h"

#include "2d/toy/gauss_kernel.h"
#include "2d/toy/reconstruction.h"
#include "2d/toy/gauss_scanner.h"

using Scanner = PET2D::Toy::GaussScanner<float>;
using Response = Scanner::Response;
using Kernel = PET2D::Toy::GaussKernel<Scanner::F>;
using Reconstrucion = PET2D::Toy::Reconstruction<Kernel, Response>;
TEST("2d/toy/reconstruction") {
  std::stringstream ss("0 0\n 0 1.1\n 1.2 0.6\n");
  Kernel kernel(1.2, 2.7);

  Reconstrucion rec;

  rec << ss;
  REQUIRE(rec.n_resp() == 3);

  REQUIRE(rec.response(1).x == Approx(0));
  REQUIRE(rec.response(1).y == Approx(1.1));
}
