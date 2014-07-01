#pragma once

#include <cmath>

#include "geometry/point.h"
#include "util/compat.h"

template <typename FType = double> class Kernel {
 public:
  typedef FType F;
  typedef FType FVec[3];
  typedef ::Point<F> Point;

  static const F INVERSE_PI = F(M_1_PI);
  static const F INVERSE_POW_TWO_PI = F(1 / (2 * M_PI * M_PI));

 private:
  static F multiply(const FVec vec_a,
                    const FVec inverse_correlation_matrix_diag,
                    const FVec vec_b) $ {
    return vec_a[0] * inverse_correlation_matrix_diag[0] * vec_b[0] +
           vec_a[1] * inverse_correlation_matrix_diag[1] * vec_b[1] +
           vec_a[2] * inverse_correlation_matrix_diag[2] * vec_b[2];
  }

 public:
  F test(const F y,
         const F z,
         const Point pixel_center,
         const F dl,
         const F sigma) $ {

    return (INVERSE_POW_TWO_PI * (1 / (sigma * dl))) *
           compat::exp(F(-0.5) *
                       (compat::pow((pixel_center.x - y) / dl, F(2)) +
                        compat::pow((pixel_center.y - z) / sigma, F(2))));
  }

  F operator()(const F y,
               const F tan,
               const F inv_cos,
               const F pow_inv_cos,
               const F R,
               const Point pixel_center,
               const FVec inverse_correlation_matrix_diag,
               const F sqrt_det_correlation_matrix) const $ {

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

    F a_ic_a = multiply(vec_a, inverse_correlation_matrix_diag, vec_a);
    F b_ic_a = multiply(vec_b, inverse_correlation_matrix_diag, vec_a);
    F b_ic_b = multiply(vec_b, inverse_correlation_matrix_diag, vec_b);
    F o_ic_b = multiply(vec_o, inverse_correlation_matrix_diag, vec_b);

    F norm = a_ic_a + (2 * o_ic_b);

    F element_before_exp =
        INVERSE_POW_TWO_PI * (sqrt_det_correlation_matrix / compat::sqrt(norm));

    F exp_element = -F(0.5) * (b_ic_b - ((b_ic_a * b_ic_a) / norm));

    F exp = compat::exp(exp_element);

    return element_before_exp * exp;
  }
};
