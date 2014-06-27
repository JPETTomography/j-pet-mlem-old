#ifndef KERNEL_H
#define KERNEL_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <utility>
#include "strip_detector.h"

#define SINGLE_INVERSE_POW_TWO_PI (1.0f / (2.0f * 3.141592f * 3.141592f))

template <typename T> class Kernel {
 public:
  static constexpr const T INVERSE_PI = T(1.0 / M_PI);
  static constexpr const T INVERSE_POW_TWO_PI = T(1.0 / (2.0 * M_PI * M_PI));
  typedef std::pair<T, T> Point;

  T multiply_elements(T* vec_a, StripDetector<T>& detector, T* vec_b) {

    T output = T(0.0);

    output += vec_a[0] * detector.inv_c(0, 0) * vec_b[0];
    output += vec_a[1] * detector.inv_c(1, 1) * vec_b[1];
    output += vec_a[2] * detector.inv_c(2, 2) * vec_b[2];

    return output;
  }

  T test_kernel(T& y, T& z, Point& pixel_center, T dl, T sigma) {

    return (SINGLE_INVERSE_POW_TWO_PI * (1 / (sigma * dl))) *
           exp(-0.5f * (pow((pixel_center.first - y) / dl, 2) +
                        pow((pixel_center.second - z) / sigma, 2)));
  }

  T calculate_kernel(T& y,
                     T& _tan,
                     T& inv_cos,
                     T& pow_inv_cos,
                     Point& pixel_center,
                     StripDetector<T>& detector,
                     T& sqrt_det_correlation_matrix) {

    T R_distance = detector.radius();
    T vec_o[3];
    T vec_a[3];
    T vec_b[3];

    vec_o[0] = -(pixel_center.first + y - R_distance) * _tan * pow_inv_cos;
    vec_o[1] = -(pixel_center.first + y + R_distance) * _tan * pow_inv_cos;
    vec_o[2] =
        -(pixel_center.first + y) * inv_cos * (T(1) + T(2) * (_tan * _tan));

    vec_a[0] = -(pixel_center.first + y - R_distance) * pow_inv_cos;
    vec_a[1] = -(pixel_center.first + y + R_distance) * pow_inv_cos;
    vec_a[2] = -T(2) * (pixel_center.first + y) * (inv_cos * _tan);

    vec_b[0] = pixel_center.second - (pixel_center.first * _tan);
    vec_b[1] = pixel_center.second - (pixel_center.first * _tan);
    vec_b[2] = -T(2) * pixel_center.first * inv_cos;

    T a_ic_a = multiply_elements(vec_a, detector, vec_a);
    T b_ic_a = multiply_elements(vec_b, detector, vec_a);
    T b_ic_b = multiply_elements(vec_b, detector, vec_b);
    T o_ic_b = multiply_elements(vec_o, detector, vec_b);

    T norm = a_ic_a + (T(2.f) * o_ic_b);

    T element_before_exp =
        INVERSE_POW_TWO_PI * (sqrt_det_correlation_matrix / std::sqrt(norm));

    T exp_element = -T(0.5f) * (b_ic_b - ((b_ic_a * b_ic_a) / norm));

    T _exp = std::exp(exp_element);

    return (element_before_exp * _exp);
  }
};

#endif  // KERNEL_H
