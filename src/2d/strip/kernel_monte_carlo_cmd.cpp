#include "2d/strip/gausian_kernel.h"

#include "common/types.h"
#include "2d/geometry/vector.h"
using Kernel = PET2D::Strip::GaussianKernel<F>;
using Vector = PET2D::Vector<F>;

int main() {

  Kernel kernel(0.01, 0.04);

  F y = 0.0;
  F R = 0.5;

  F d = 0.01;
  F dfi = 0.01;

  double sum = 0.0;
  for (F angle = -M_PI / 4; angle < M_PI / 4; angle += dfi) {
    F tan = std::tan(angle);
    F sec = 1 / cos(angle);
    for (F dx = -0.2; dx <= 0.2; dx += d) {
      for (F dy = -0.2; dy <= 0.2; dy += d) {
        sum += kernel(y + dy, tan, sec, R, Vector(dx, dy));
      }
    }
  }
  std::cout << sum * d * d * dfi << "\n";
  return 0;
}
