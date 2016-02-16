#ifndef GAUSS_KERNEL
#define GAUSS_KERNEL

#include <cmath>

namespace PET1D {
namespace Toy {

template <typename FType> class GaussKernel {
 public:
  using F = FType;

  GaussKernel(F sigma)
      : sigma(sigma),
        norm(1 / (std::sqrt(2 * M_PI) * sigma)),
        inv_two_sigma2_(1 / (2 * sigma * sigma)) {}

  F operator()(F e, F x) {
    F diff = e - x;
    return norm * std::exp(-diff * diff * inv_two_sigma2_);
  }

  F sigma;
  F norm;

 private:
  F inv_two_sigma2_;
};
}
}
#endif  // GAUSS_KERNEL
