#include "2d/strip/gausian_kernel.h"

#include "common/types.h"
#include "2d/geometry/vector.h"
using Kernel = PET2D::Strip::GaussianKernel<F>;
using Vector = PET2D::Vector<F>;

int main() {

  Kernel kernel(0.01, 0.04);

  F y = 0.1;
  F R = 0.5;
  F angle= M_PI/4;
  F tan = std::tan(angle);
  F sec = 1/cos(angle);

  F d = 0.01;
  double sum = 0.0;
  for (F rx = -0.2; rx <= 0.2; rx += d) {
    for (F ry = -0.2; ry <= 0.2; ry += d) {
      sum += kernel(y, tan, sec, R, Vector(rx, ry));
    }
  }

  std::cout << sum * d * d << "\n";
  return 0;
}
