#pragma once

#include <cmath>

#include "strip_detector.h"

template <typename FType = double> class Kernel {
 public:
  typedef FType F;
  static constexpr const F INVERSE_PI = F(M_1_PI);
  static constexpr const F INVERSE_POW_TWO_PI = F(1 / (2 * M_PI * M_PI));
  typedef std::pair<F, F> Point;

  F multiply_elements(F* vec_a, StripDetector<F>& detector, F* vec_b) {

    F output = 0;

    output += vec_a[0] * detector.inv_c(0, 0) * vec_b[0];
    output += vec_a[1] * detector.inv_c(1, 1) * vec_b[1];
    output += vec_a[2] * detector.inv_c(2, 2) * vec_b[2];

    return output;
  }

  F test_kernel(F& y, F& z, Point& pixel_center, F dl, F sigma) {

    return (INVERSE_POW_TWO_PI * (1 / (sigma * dl))) *
           std::exp(F(-0.5) * (std::pow((pixel_center.first - y) / dl, 2) +
                               std::pow((pixel_center.second - z) / sigma, 2)));
  }

  F calculate_kernel(F& y,
                     F& tan,
                     F& inv_cos,
                     F& pow_inv_cos,
                     Point& pixel_center,
                     StripDetector<F>& detector,
                     F& sqrt_det_correlation_matrix) {

    F R_distance = detector.radius;
    F vec_o[3];
    F vec_a[3];
    F vec_b[3];

    vec_o[0] = -(pixel_center.first + y - R_distance) * tan * pow_inv_cos;
    vec_o[1] = -(pixel_center.first + y + R_distance) * tan * pow_inv_cos;
    vec_o[2] = -(pixel_center.first + y) * inv_cos * (1 + 2 * (tan * tan));

    vec_a[0] = -(pixel_center.first + y - R_distance) * pow_inv_cos;
    vec_a[1] = -(pixel_center.first + y + R_distance) * pow_inv_cos;
    vec_a[2] = -2 * (pixel_center.first + y) * (inv_cos * tan);

    vec_b[0] = pixel_center.second - (pixel_center.first * tan);
    vec_b[1] = pixel_center.second - (pixel_center.first * tan);
    vec_b[2] = -2 * pixel_center.first * inv_cos;

    F a_ic_a = multiply_elements(vec_a, detector, vec_a);
    F b_ic_a = multiply_elements(vec_b, detector, vec_a);
    F b_ic_b = multiply_elements(vec_b, detector, vec_b);
    F o_ic_b = multiply_elements(vec_o, detector, vec_b);

    F norm = a_ic_a + (2 * o_ic_b);

    F element_before_exp =
        INVERSE_POW_TWO_PI * (sqrt_det_correlation_matrix / std::sqrt(norm));

    F exp_element = -F(0.5) * (b_ic_b - ((b_ic_a * b_ic_a) / norm));

    F _exp = std::exp(exp_element);

    return (element_before_exp * _exp);
  }
};
