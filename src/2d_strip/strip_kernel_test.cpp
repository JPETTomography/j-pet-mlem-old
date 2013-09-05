#include "catch.hpp"
#include <vector>
#include <cmath>

typedef std::pair<int, int> Pixel;
typedef std::pair<float, float> Point;

template <typename T>
T multiply_elements(std::vector<T>& vec_a,
                    std::vector<T>& vec_b,
                    std::vector<std::vector<T>>& inverse_correlation_matrix) {

  std::vector<T> a(vec_a.size(), T(0));

  float output = 0.f;
  // add AVX
  for (unsigned i = 0; i < vec_a.size(); ++i) {
    for (unsigned j = 0; j < vec_b.size(); j++) {
      a[i] += vec_a[j] * inverse_correlation_matrix[j][i];
    }
    output += a[i] * vec_b[i];
  }

  return output;
}

template <typename T> T kernel(T y, T angle, Point pixel_center) {

  T pow_sigma_z = 10 * 10;
  T pow_sigma_dl = 63 * 63;
  T R_distance = 500;

  std::vector<std::vector<T>> inverse_correlation_matrix;

  inverse_correlation_matrix.resize(3, std::vector<T>(3, T()));

  inverse_correlation_matrix[0][0] = T(1) / pow_sigma_z;
  inverse_correlation_matrix[0][1] = T(0.0f);
  inverse_correlation_matrix[0][2] = T(0.0f);

  inverse_correlation_matrix[1][0] = T(0.0f);
  inverse_correlation_matrix[1][1] = T(1) / pow_sigma_z;
  inverse_correlation_matrix[1][2] = T(0.0f);

  inverse_correlation_matrix[2][0] = T(0.0f);
  inverse_correlation_matrix[2][1] = T(0.0f);
  inverse_correlation_matrix[2][2] = T(1) / pow_sigma_dl;

  T sqrt_det_correlation_matrix = std::sqrt(inverse_correlation_matrix[0][0] *
                                            inverse_correlation_matrix[1][1] *
                                            inverse_correlation_matrix[2][2]);

  T _tan = std::tan(angle);
  T inv_cos = T(1) / std::cos(angle);
  T pow_inv_cos = inv_cos * inv_cos;

  std::vector<T> vec_o(3, T());
  std::vector<T> vec_a(3, T());
  std::vector<T> vec_b(3, T());

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

  T a_ic_a = multiply_elements(vec_a, vec_a, inverse_correlation_matrix);
  T b_ic_a = multiply_elements(vec_b, vec_a, inverse_correlation_matrix);
  T b_ic_b = multiply_elements(vec_b, vec_b, inverse_correlation_matrix);
  T o_ic_b = multiply_elements(vec_o, vec_b, inverse_correlation_matrix);

  T norm = a_ic_a + (T(2) * o_ic_b);

  T element_before_exp = (1.0 / (2.0 * M_PI * M_PI)) *
                         (sqrt_det_correlation_matrix / std::sqrt(norm));

  T exp_element = -T(0.5) * (b_ic_b - ((b_ic_a * b_ic_a) / norm));

  T _exp = std::exp(exp_element);

  return (element_before_exp * _exp);
}

TEST_CASE("kernel tests", "kernel") {

  // POINT = -25.35f,25.35f
  CHECK(kernel<float>(-25.35f, 0.65387f, Point(-22.5f, 27.5f)) ==
        Approx(1.15581e-16));

  // CHECK(vec_o[1] == Approx(-2100.0));
  // CHECK(vec_o[2] == Approx(-1800 * std::sqrt(2)));
}
