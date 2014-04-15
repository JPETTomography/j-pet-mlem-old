#include <cuda_runtime.h>

#define SINGLE_INVERSE_PI 0.318309f
#define SINGLE_INVERSE_POW_TWO_PI (1.0f / (2.0f * 3.141592f * 3.141592f))

/*
float inv_c(int i, int j) const { return inverse_correlation_matrix_[i][j]; }

event<F> to_projection_space_tan(const ImageSpaceEventTan<F> &ev) {
  F z_u = ev.z + (radius() - ev.y) * ev.tan;
  F z_d = ev.z - (radius() + ev.y) * ev.tan;
  F dl = -F(2.0) * ev.y * sqrt(ev.tan * ev.tan + F(1.0));
  return event<F>(z_u, z_d, dl);
}



Point pixel_center(int i, int j) {
  return std::make_pair<F>(
      grid_ul_y_ - i * pixel_height_ - 0.5 * pixel_height_,
      grid_ul_z_ + j * pixel_width_ + 0.5 * pixel_width_);
}

Point pixel_center(Pixel pix) { return pixel_center(pix.first, pix.second); }
*/

__device__ inline float multiply_elements(float* vec_a,
                                          volatile float* inv_c,
                                          float* vec_b) {

  float output = 0.0f;

  output += vec_a[0] * inv_c[0] * vec_b[0];
  output += vec_a[1] * inv_c[1] * vec_b[1];
  output += vec_a[2] * inv_c[2] * vec_b[2];

  return output;
}

__device__ float calculate_kernel(float& y,
                                  float& _tan,
                                  float& inv_cos,
                                  float& pow_inv_cos,
                                  float2& pixel_center,
                                  float* inv_c,
                                  gpu_config::GPU_parameters& cfg,
                                  float& sqrt_det_correlation_matrix) {

  float vec_o[3];
  float vec_a[3];
  float vec_b[3];

  vec_o[0] = -(pixel_center.x + y - cfg.R_distance) * _tan * pow_inv_cos;
  vec_o[1] = -(pixel_center.x + y + cfg.R_distance) * _tan * pow_inv_cos;
  vec_o[2] = -(pixel_center.x + y) * inv_cos * (1.0f + 2.0f * (_tan * _tan));

  vec_a[0] = -(pixel_center.x + y - cfg.R_distance) * pow_inv_cos;
  vec_a[1] = -(pixel_center.x + y + cfg.R_distance) * pow_inv_cos;
  vec_a[2] = -2.0f * (pixel_center.x + y) * (inv_cos * _tan);

  vec_b[0] = pixel_center.y - (pixel_center.x * _tan);
  vec_b[1] = pixel_center.y - (pixel_center.x * _tan);
  vec_b[2] = -2.0 * pixel_center.x * inv_cos;

  float a_ic_a = multiply_elements(vec_a, inv_c, vec_a);
  float b_ic_a = multiply_elements(vec_b, inv_c, vec_a);
  float b_ic_b = multiply_elements(vec_b, inv_c, vec_b);
  float o_ic_b = multiply_elements(vec_o, inv_c, vec_b);

  float norm = a_ic_a + (2.f * o_ic_b);

#if TEST
  if (threadIdx.x == 0 && blockIdx.x == 0) {

    printf("vec_o[0]: %f\n", vec_o[0]);
    printf("vec_o[1]: %f\n", vec_o[1]);
    printf("vec_o[2]: %f\n", vec_o[2]);

    printf("vec_a[0]: %f\n", vec_a[0]);
    printf("vec_a[1]: %f\n", vec_a[1]);
    printf("vec_a[2]: %f\n", vec_a[2]);

    printf("vec_b[0]: %f\n", vec_b[0]);
    printf("vec_b[1]: %f\n", vec_b[1]);
    printf("vec_b[2]: %f\n", vec_b[2]);

    printf("inv_c[0]: %f\n", inv_c[0]);
    printf("inv_c[1]: %f\n", inv_c[1]);
    printf("inv_c[2]: %f\n", inv_c[2]);

    printf("a_ic_a: %f\n", a_ic_a);
    printf("b_ic_a: %f\n", b_ic_a);
    printf("b_ic_b: %f\n", b_ic_b);
    printf("o_ic_b: %f\n", o_ic_b);
    printf("norm: %f\n", norm);

    float element_before_exp = SINGLE_INVERSE_POW_TWO_PI *
                               (sqrt_det_correlation_matrix / std::sqrt(norm));

    printf("/------------------------------------------/\n");

    printf("SINGLE_INVERSE_POW_TWO_PI: %f\n", SINGLE_INVERSE_POW_TWO_PI);
    printf("sqrt_det_correlation_matrix: %f\n", sqrt_det_correlation_matrix);

    printf("/------------------------------------------/\n");

    printf("Element before exp: %ef\n", element_before_exp);
    printf("Exp: %ef \n", exp(-(0.5f) * (b_ic_b - ((b_ic_a * b_ic_a) / norm))));
  }

#endif
  return (SINGLE_INVERSE_POW_TWO_PI *
          (sqrt_det_correlation_matrix / sqrt(norm)) *
          exp(-(0.5f) * (b_ic_b - ((b_ic_a * b_ic_a) / norm))));
}

__device__ float2 pixel_center(int& i,
                               int& j,
                               float& pixel_height_,
                               float& pixel_width_,
                               float& grid_size_y_,
                               float& grid_size_z_) {

  return make_float2(
      0.5f * grid_size_y_ - i * pixel_height_ - 0.5f * pixel_height_,
      -0.5f * grid_size_z_ + j * pixel_width_ + 0.5f * pixel_width_);
}

__device__ float2 pixel_location(float& y,
                                 float& z,
                                 float& pixel_height_,
                                 float& pixel_width_,
                                 float& grid_size_y_,
                                 float& grid_size_z_) {

  return make_float2(floor((0.5f * grid_size_y_) - y / pixel_height_),
                     floor(z + (0.5 * grid_size_z_) / pixel_width_));
}

// Pixel pixel_location(float2& p) { return pixel_location(p.x, p.y); }

__device__ float sensitivity(float& y,
                             float& z,
                             float& radius,
                             float& half_scintilator_length) {
  float L_plus = (half_scintilator_length + z);
  float L_minus = (half_scintilator_length - z);
  float R_plus = radius + y;
  float R_minus = radius - y;
  return SINGLE_INVERSE_PI *
         (std::atan(min(L_minus / R_minus, L_plus / R_plus)) -
          std::atan(max(-L_plus / R_minus, -L_minus / R_plus)));
}

__device__ float bbz(float& A, float& C, float& B_2) {
  return (3.0f) / sqrt(C - (B_2 / A));
}

__device__ float bby(float& A, float& C, float& B_2) {
  return (3.0f) / sqrt(A - (B_2 / C));
}

__device__ int pixels_in_line(float& length, float& pixel_size) {
  float d = (length + 0.5f) / pixel_size;
  return int(d);
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
