#pragma once

#include <cmath>
#if _MSC_VER
#include <array>
#endif

#include "geometry/point.h"
#include "util/compat.h"

template <typename FType = double> class Kernel {
 public:
  typedef FType F;
  typedef ::Point<F> Point;
#if !_MSC_VER || __CUDACC__
  typedef FType FVec[3];
#else
  typedef std::array<F, 3> FVec;
#endif

 private:
  _ static F multiply(const FVec vec_a,
                      const FVec inv_cor_mat_diag,
                      const FVec vec_b) {
    return vec_a[0] * inv_cor_mat_diag[0] * vec_b[0] +
           vec_a[1] * inv_cor_mat_diag[1] * vec_b[1] +
           vec_a[2] * inv_cor_mat_diag[2] * vec_b[2];
  }

 public:
  _ F test(const F y,
           const F z,
           const Point pixel_center,
           const F dl,
           const F sigma) {

    const F INV_POW_TWO_PI = F(1 / (2 * M_PI * M_PI));

    return (INV_POW_TWO_PI * (1 / (sigma * dl))) *
           compat::exp(F(-0.5) *
                       (compat::pow((pixel_center.x - y) / dl, F(2)) +
                        compat::pow((pixel_center.y - z) / sigma, F(2))));
  }

  _ F operator()(const F y,
                 const F tan,
                 const F inv_cos,
                 const F pow_inv_cos,
                 const F R,
                 const Point pixel_center,
                 const FVec inv_cor_mat_diag,
                 const F sqrt_det_cor_mat) const {

    FVec vec_o;
    FVec vec_a;
    FVec vec_b;

    vec_o[0] = -(pixel_center.x + y - R) * tan * pow_inv_cos;
    vec_o[1] = -(pixel_center.x + y + R) * tan * pow_inv_cos;
    vec_o[2] = -(pixel_center.x + y) * inv_cos * (1 + 2 * (tan * tan));

    vec_a[0] = -(pixel_center.x + y - R) * pow_inv_cos;
    vec_a[1] = -(pixel_center.x + y + R) * pow_inv_cos;
    vec_a[2] = -2 * (pixel_center.x + y) * (inv_cos * tan);

    vec_b[0] = pixel_center.y - (pixel_center.x * tan);
    vec_b[1] = pixel_center.y - (pixel_center.x * tan);
    vec_b[2] = -2 * pixel_center.x * inv_cos;

    F a_ic_a = multiply(vec_a, inv_cor_mat_diag, vec_a);
    F b_ic_a = multiply(vec_b, inv_cor_mat_diag, vec_a);
    F b_ic_b = multiply(vec_b, inv_cor_mat_diag, vec_b);
    F o_ic_b = multiply(vec_o, inv_cor_mat_diag, vec_b);

    F norm = a_ic_a + (2 * o_ic_b);

    const F INV_POW_TWO_PI = F(1 / (2 * M_PI * M_PI));

    F element_before_exp =
        INV_POW_TWO_PI * (sqrt_det_cor_mat / compat::sqrt(norm));

    F exp_element = -F(0.5) * (b_ic_b - ((b_ic_a * b_ic_a) / norm));

    F exp = compat::exp(exp_element);

    return element_before_exp * exp;
  }
};
