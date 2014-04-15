#pragma once

#include <cuda_runtime.h>
#include "../config.h"
#include "../event.h"
#include "reconstruction_methods.cuh"

__device__ float lookup_table[200][200];

__shared__ float shmem_z_u[128];
__shared__ float shmem_z_d[128];
__shared__ float shmem_delta_l[128];

template <typename T>
__global__ void reconstruction_2d_strip_cuda(gpu_config::GPU_parameters cfg,
                                             event<float>* event_list,
                                             int event_list_size,
                                             float* image_buffor,
                                             float* rho) {

  int tid = (blockIdx.x * blockDim.x) + threadIdx.x;

  // future constant memory allocation on host side !FIXME

  __shared__ float inv_c[3];

  if (threadIdx.x == 0) {

    inv_c[0] = cfg.inv_pow_sigma_z;
    inv_c[1] = cfg.inv_pow_sigma_z;
    inv_c[2] = cfg.inv_pow_sigma_dl;
  }

  for (int i = 0; i < event_list_size / blockDim.x; ++i) {
    float z_u = event_list[tid + cfg.number_of_blocks * blockDim.x].z_u;
    float z_d = event_list[tid + cfg.number_of_blocks * blockDim.x].z_d;
    float delta_l = event_list[tid + cfg.number_of_blocks * blockDim.x].dl;
    float acc = 0.f;

    float sqrt_det_correlation_matrix =
        sqrt(cfg.inv_pow_sigma_z * cfg.inv_pow_sigma_z * cfg.inv_pow_sigma_dl);

    // angle space transformation
    float tn = event_tan(z_u, z_d, cfg.R_distance);
    float y = event_y(delta_l, tn);
    float z = event_z(z_u, z_d, y, tn);

    float cos_ = cos(tn);

    float sec_ = float(1.0f) / cos_;
    float sec_sq_ = sec_ * sec_;

    float A = (((T(4.0f) / (cos_ * cos_)) * cfg.inv_pow_sigma_dl) +
               (T(2.0f) * tn * tn * cfg.inv_pow_sigma_z));
    float B = -T(4.0f) * tn * cfg.inv_pow_sigma_z;
    float C = T(2.0f) * cfg.inv_pow_sigma_z;
    float B_2 = (B / T(2.0f)) * (B / T(2.0f));

    float bb_y = bby(A, C, B_2);

    float bb_z = bbz(A, C, B_2);

    float2 center_pixel = make_float2(y, z);

#if BOUNDING_BOX_TEST

    // check(0.0470448, 94.3954, 21.6737);
    // check(-0.594145, 78.3053, 56.9959);
    // check(0.20029, 92.6108, 28.3458);
    // check(-0.571667, 79.4745, 55.3539);
    // check(-0.420542, 86.266, 44.0276);

    if (tid == 0) {

      float angle = -0.594145f;

      float test_A =
          (((4.0f / (std::cos(angle) * std::cos(angle))) *
            cfg.inv_pow_sigma_dl) +
           (2.0f * std::tan(angle) * std::tan(angle) * cfg.inv_pow_sigma_z));
      float test_B = -4.0f * std::tan(angle) * cfg.inv_pow_sigma_z;
      float test_C = 2.0f * cfg.inv_pow_sigma_z;
      float test_B_2 = (test_B / 2.0f) * (test_B / 2.0f);

      float test_bb_y = bby(test_A, test_C, test_B_2);

      float test_bb_z = bbz(test_A, test_C, test_B_2);

      printf(
          "sig_dl: %f sig_z: %f\n", cfg.inv_pow_sigma_dl, cfg.inv_pow_sigma_z);
      printf("BB: %f %f \n", test_bb_y, test_bb_z);
    }
#endif

    // bounding box limits for event
    float2 ur =
        make_float2(center_pixel.x - pixels_in_line(bb_y, cfg.pixel_size),
                    center_pixel.y + pixels_in_line(bb_z, cfg.pixel_size));
    float2 dl =
        make_float2(center_pixel.x + pixels_in_line(bb_y, cfg.pixel_size),
                    center_pixel.y - pixels_in_line(bb_z, cfg.pixel_size));
    float2 pp;

    for (int iz = dl.y; iz < ur.y; ++iz) {
      for (int iy = ur.x; iy < dl.x; ++iy) {

        pp = pixel_center(iy,
                          iz,
                          cfg.pixel_size,
                          cfg.pixel_size,
                          cfg.grid_size_y_,
                          cfg.grid_size_z_);

        if (in_ellipse(A, B, C, y, z, pp)) {

          pp.x -= y;
          pp.y -= z;

          T event_kernel = calculate_kernel(y,
                                            tn,
                                            sec_,
                                            sec_sq_,
                                            pp,
                                            inv_c,
                                            cfg,
                                            sqrt_det_correlation_matrix) /
                           lookup_table[iy][iz];

          acc += event_kernel * lookup_table[iy][iz] *
                 rho[iy + (iz * cfg.n_pixels)];
        }
      }
    }

    for (int iz = dl.y; iz < ur.y; ++iz) {
      for (int iy = ur.x; iy < dl.x; ++iy) {

        pp = pixel_center(iy,
                          iz,
                          cfg.pixel_size,
                          cfg.pixel_size,
                          cfg.grid_size_y_,
                          cfg.grid_size_z_);

        if (in_ellipse(A, B, C, y, z, pp)) {

          pp.x -= y;
          pp.y -= z;

          T event_kernel = calculate_kernel(y,
                                            tn,
                                            sec_,
                                            sec_sq_,
                                            pp,
                                            inv_c,
                                            cfg,
                                            sqrt_det_correlation_matrix) /
                           lookup_table[iy][iz];

          atomicExch(&image_buffor[(blockIdx.x * cfg.n_pixels * cfg.n_pixels) +
                                   iy + (iz * cfg.n_pixels)],
                     image_buffor[(blockIdx.x * cfg.n_pixels * cfg.n_pixels) +
                                  iy + (iz * cfg.n_pixels)] +
                         (event_kernel * rho[iy + (iz * cfg.n_pixels)] / acc));
        }
      }
    }
  }
}
