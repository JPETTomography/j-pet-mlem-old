#pragma once

#include <cuda_runtime.h>
#include <cmath>

#include "geometry/pixel.h"
#include "geometry/point.h"

#include "config.h"

static const float INVERSE_PI = (float)M_1_PI;
static const float INVERSE_POW_TWO_PI = (float)(1 / (2 * M_PI * M_PI));

#if OLD_WARP_SPACE_PIXEL
__device__ void warp_space_pixel(Pixel<>& pixel,
                                 int offset,
                                 int ul,
                                 int width,
                                 int height,
                                 int& index) {

  index = threadIdx.x & 31 + offset;
  pixel.y = index / width;
  pixel.x = index - width * pixel.y;
  pixel.y += ul.y;
  pixel.x += ul.x;
}
#else
__device__ Pixel<> warp_space_pixel(int& offset,
                                    Pixel<>& first_pixel,
                                    Pixel<> ul,
                                    Pixel<> ur,
                                    Pixel<> dl,
                                    int tid) {
  offset = ur.y - first_pixel.y;
  int tid_in_warp = threadIdx.x & 31;

  Pixel<> pixel;
  if (tid_in_warp < offset) {
    pixel = Pixel<>(first_pixel.x, first_pixel.y + tid_in_warp);
  } else {
    pixel = Pixel<>(
        first_pixel.x + (int)(__fdividef(tid - offset, ur.y - dl.y)) + 1,
        dl.y + ((tid - offset) % (ur.y - dl.y)));
  }

  first_pixel.y = ul.y + (32 - offset) % (ur.y - ul.y);
  first_pixel.x += int(__fdividef(32 - offset, ur.y - ul.y)) + 1;

  return pixel;
}
#endif

template <typename F>
__device__ float multiply(F* vec_a, volatile F* inv_c, F* vec_b) {
  return vec_a[0] * inv_c[0] * vec_b[0] +  //
         vec_a[1] * inv_c[1] * vec_b[1] +  //
         vec_a[2] * inv_c[2] * vec_b[2];
}

template <typename F>
__device__ float main_kernel(F y,
                             F tan,
                             F inv_cos,
                             F pow_inv_cos,
                             Point<F> pixel_center,
                             F* inv_c,
                             CUDA::Config<F>& cfg,
                             F sqrt_det_correlation_matrix) {

  F vec_o[3];
  F vec_a[3];
  F vec_b[3];

  vec_o[0] = -(pixel_center.x + y - cfg.R_distance) * tan * pow_inv_cos;
  vec_o[1] = -(pixel_center.x + y + cfg.R_distance) * tan * pow_inv_cos;
  vec_o[2] = -(pixel_center.x + y) * inv_cos * (1 + 2 * (tan * tan));

  vec_a[0] = -(pixel_center.x + y - cfg.R_distance) * pow_inv_cos;
  vec_a[1] = -(pixel_center.x + y + cfg.R_distance) * pow_inv_cos;
  vec_a[2] = -2 * (pixel_center.x + y) * (inv_cos * tan);

  vec_b[0] = pixel_center.y - (pixel_center.x * tan);
  vec_b[1] = pixel_center.y - (pixel_center.x * tan);
  vec_b[2] = -2 * pixel_center.x * inv_cos;

  F a_ic_a = multiply<F>(vec_a, inv_c, vec_a);
  F b_ic_a = multiply<F>(vec_b, inv_c, vec_a);
  F b_ic_b = multiply<F>(vec_b, inv_c, vec_b);
  F o_ic_b = multiply<F>(vec_o, inv_c, vec_b);

  F norm = a_ic_a + (2.f * o_ic_b);

  return INVERSE_POW_TWO_PI * (sqrt_det_correlation_matrix / sqrt(norm)) *
         exp(F(-0.5) * (b_ic_b - ((b_ic_a * b_ic_a) / norm)));
}

template <typename F>
__device__ float test_kernel(F y,
                             F z,
                             Point<F> pixel_center,
                             CUDA::Config<F>& cfg) {

  return INVERSE_POW_TWO_PI * (1 / (cfg.sigma * cfg.dl)) *
         exp(F(-0.5) * (pow((pixel_center.x - y) / cfg.dl, 2) +
                        pow((pixel_center.y - z) / cfg.sigma, 2)));
}

template <typename F>
__host__ __device__ Point<F> pixel_center(int y,
                                          int z,
                                          F pixel_height,
                                          F pixel_width,
                                          F half_grid_size,
                                          F half_pixel_size) {

  return Point<F>(half_grid_size - y * pixel_height - half_pixel_size,
                  -half_grid_size + z * pixel_width + half_pixel_size);
}

__device__ Pixel<> pixel_location(float y,
                                  float z,
                                  float pixel_height,
                                  float pixel_width,
                                  float grid_size_y,
                                  float grid_size_z) {

  return Pixel<>(floor(((0.5f * grid_size_y) - y) / pixel_height),
                 floor((z - (-0.5f * grid_size_z)) / pixel_width));
}

__host__ __device__ float sensitivity(float y,
                                      float z,
                                      float radius,
                                      float half_scintilator_length) {

  float L_plus = (half_scintilator_length + z);
  float L_minus = (half_scintilator_length - z);
  float R_plus = radius + y;
  float R_minus = radius - y;

  return INVERSE_PI * (atanf(min(L_minus / R_minus, L_plus / R_plus)) -
                       atanf(max(-L_plus / R_minus, -L_minus / R_plus)));
}

__device__ float bbz(float A, float C, float B_2) {
  return 3 / sqrtf(C - (B_2 / A));
}

__device__ float bby(float A, float C, float B_2) {
  return 3 / sqrtf(A - (B_2 / C));
}

__device__ int pixels_in_line(float length, float pixel_size) {
  return int((length + 0.5f) / pixel_size);
}

template <typename F>
__device__ bool in_ellipse(F A, F B, F C, F y, F z, Point<F> point) {

  F dy = (point.x - y);
  F dz = (point.y - z);

  // quadratic ellipse equation check
  return (A * (dy * dy)) + (B * dy * dz) + (C * (dz * dz)) <= 9;
}
