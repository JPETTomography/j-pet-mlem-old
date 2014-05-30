#include <cuda_runtime.h>

#define SINGLE_INVERSE_PI 0.3183098861f
#define SINGLE_INVERSE_POW_TWO_PI (1.0f / (2.0f * 3.141592f * 3.141592f))

#define DEBUG_KERNEL 0

template <typename T>
__device__ float multiply_elements(T* vec_a, volatile T* inv_c, T* vec_b) {

  T output = 0.0f;

  output +=  vec_a[0] * inv_c[0] * vec_b[0];
  output += vec_a[1] * inv_c[1] * vec_b[1];
  output += vec_a[2] * inv_c[2] * vec_b[2];

  return output;
}

template <typename T>
__device__ float calculate_kernel(T& y,
                                  T& _tan,
                                  T& inv_cos,
                                  T& pow_inv_cos,
                                  float2& pixel_center,
                                  T* inv_c,
                                  gpu_config::GPU_parameters& cfg,
                                  T& sqrt_det_correlation_matrix) {

  T vec_o[3];
  T vec_a[3];
  T vec_b[3];

  vec_o[0] = -(pixel_center.x + y - cfg.R_distance) * _tan * pow_inv_cos;
  vec_o[1] = -(pixel_center.x + y + cfg.R_distance) * _tan * pow_inv_cos;
  vec_o[2] = -(pixel_center.x + y) * inv_cos * (1.0f + 2.0f * (_tan * _tan));

  vec_a[0] = -(pixel_center.x + y - cfg.R_distance) * pow_inv_cos;
  vec_a[1] = -(pixel_center.x + y + cfg.R_distance) * pow_inv_cos;
  vec_a[2] = -2.0f * (pixel_center.x + y) * (inv_cos * _tan);

  vec_b[0] = pixel_center.y - (pixel_center.x * _tan);
  vec_b[1] = pixel_center.y - (pixel_center.x * _tan);
  vec_b[2] = -2.0f * pixel_center.x * inv_cos;

  T a_ic_a = multiply_elements<T>(vec_a, inv_c, vec_a);
  T b_ic_a = multiply_elements<T>(vec_b, inv_c, vec_a);
  T b_ic_b = multiply_elements<T>(vec_b, inv_c, vec_b);
  T o_ic_b = multiply_elements<T>(vec_o, inv_c, vec_b);

  T norm = a_ic_a + (2.f * o_ic_b);

  return (SINGLE_INVERSE_POW_TWO_PI *
          (sqrt_det_correlation_matrix / sqrt(norm)) *
          exp(-(0.5f) * (b_ic_b - ((b_ic_a * b_ic_a) / norm))));
}

__host__ __device__ float2 pixel_center(int& i,
                                        int& j,
                                        float& pixel_height_,
                                        float& pixel_width_,
                                        float& grid_size_y_,
                                        float& grid_size_z_,
                                        float& half_grid_size,
                                        float& half_pixel_size) {

  return make_float2(half_grid_size - i * pixel_height_ - half_pixel_size,
                     -half_grid_size + j * pixel_width_ + half_pixel_size);
}

__device__ int2 pixel_location(float& y,
                               float& z,
                               float& pixel_height_,
                               float& pixel_width_,
                               float& grid_size_y_,
                               float& grid_size_z_) {

  return make_int2(floor(((0.5f * grid_size_y_) - y) / pixel_height_),
                   floor((z - (-0.5f * grid_size_z_)) / pixel_width_));
}

__host__ __device__ float sensitivity(float y,
                                      float z,
                                      float radius,
                                      float half_scintilator_length) {

  float L_plus = (half_scintilator_length + z);
  float L_minus = (half_scintilator_length - z);
  float R_plus = radius + y;
  float R_minus = radius - y;

  return SINGLE_INVERSE_PI * (atanf(min(L_minus / R_minus, L_plus / R_plus)) -
                              atanf(max(-L_plus / R_minus, -L_minus / R_plus)));
}

__device__ float bbz(float& A, float& C, float& B_2) {
  return (3.0f) / sqrtf(C - (B_2 / A));
}

__device__ float bby(float& A, float& C, float& B_2) {
  return (3.0f) / sqrtf(A - (B_2 / C));
}

__device__ int pixels_in_line(float& length, float& pixel_size) {
  return int((length + 0.5f) / pixel_size);
}

__device__ float event_tan(float& z_u, float& z_d, float& R) {
  return (z_u - z_d) / ((2.0f) * R);
}
__device__ float event_y(float& dl, float& tan_event) {
  return -(0.5f) * (dl / sqrt((1.0f) + (tan_event * tan_event)));
}
__device__ float event_z(float& z_u, float& z_d, float& y, float& tan_event) {
  return (0.5f) * (z_u + z_d + ((2.0f) * y * tan_event));
}

__device__ bool in_ellipse(float& A,
                           float& B,
                           float& C,
                           float& y,
                           float& z,
                           float2 point) {

  float dy = (point.x - y);
  float dz = (point.y - z);

  // quadratic ellipse equation check
  return (((A * (dy * dy)) + (B * dy * dz) + (C * (dz * dz)))) <= (9.0f)
             ? true
             : false;
}