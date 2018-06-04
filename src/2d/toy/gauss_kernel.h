#pragma once

/// Two-dimensional toy (experiment) space
namespace PET2D {
/// Simple PET working in two-dimensional space
namespace Toy {

template <typename FType> class GaussKernel {

 public:
  using F = FType;

  GaussKernel(F sigma_x, F sigma_y) : sigma_x(sigma_x), sigma_y(sigma_y) {
    norm_ = 1.0 / (2 * M_PI * sigma_x * sigma_y);
    inv_2_sigma2_x_ = 1.0 / (2 * sigma_x * sigma_x);
    inv_2_sigma2_y_ = 1.0 / (2 * sigma_y * sigma_y);
  }

  F operator()(F x, F y) {
    return norm_ * exp(-inv_2_sigma2_x_ * x * x - inv_2_sigma2_y_ * y * y);
  }

  const F sigma_x;
  const F sigma_y;

  F three_sigma_bb_x() const { return 3.0 * sigma_x; }

  F three_sigma_bb_y() const { return 3.0 * sigma_y; }

 private:
  F inv_2_sigma2_x_;
  F inv_2_sigma2_y_;
  F norm_;
};
}
}
