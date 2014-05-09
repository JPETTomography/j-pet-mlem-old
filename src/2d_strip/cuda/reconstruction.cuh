#pragma once

#include <cuda_runtime.h>
#include "../config.h"
#include "../event.h"
#include "reconstruction_methods.cuh"

#define IMAGE_SPACE_LINEAR_INDEX(Y, Z) (Y * cfg.n_pixels) + Z
#define BUFFOR_LINEAR_INDEX(Y, Z) \
  (blockIdx.x * cfg.n_pixels * cfg.n_pixels) + (Y * cfg.n_pixels) + Z

template <typename T>
__global__ void reconstruction_2d_strip_cuda(gpu_config::GPU_parameters cfg,
                                             event<T>* event_list,
                                             int event_list_size,
                                             float* image_buffor,
                                             float* rho,
                                             cudaTextureObject_t tex) {

  int tid = (blockIdx.x * blockDim.x) + threadIdx.x;

  __shared__ float inv_c[3];

  int block_chunk =
      int(ceilf(event_list_size /
                (cfg.number_of_blocks * cfg.number_of_threads_per_block)));

  if (threadIdx.x == 0) {

    inv_c[0] = cfg.inv_pow_sigma_z;
    inv_c[1] = cfg.inv_pow_sigma_z;
    inv_c[2] = cfg.inv_pow_sigma_dl;
  }

  __syncthreads();

  for (int i = 0; i < block_chunk; ++i) {

    if ((i * cfg.number_of_blocks * cfg.number_of_threads_per_block) + tid <
        event_list_size) {

      float z_u = event_list[(i * cfg.number_of_blocks *
                              cfg.number_of_threads_per_block) +
                             tid].z_u;
      float z_d = event_list[(i * cfg.number_of_blocks *
                              cfg.number_of_threads_per_block) +
                             tid].z_d;
      float delta_l = event_list[(i * cfg.number_of_blocks *
                                  cfg.number_of_threads_per_block) +
                                 tid].dl;
      float acc = 0.f;

      float sqrt_det_correlation_matrix = sqrt(
          cfg.inv_pow_sigma_z * cfg.inv_pow_sigma_z * cfg.inv_pow_sigma_dl);

      // angle space transformation
      float tn = event_tan(z_u, z_d, cfg.R_distance);
      float y = event_y(delta_l, tn);
      float z = event_z(z_u, z_d, y, tn);

      float angle = atan(tn);

      float cos_ = __cosf(angle);

      float sec_ = float(1.0f) / cos_;
      float sec_sq_ = sec_ * sec_;

      float A = (((T(4.0f) / (cos_ * cos_)) * cfg.inv_pow_sigma_dl) +
                 (T(2.0f) * tn * tn * cfg.inv_pow_sigma_z));
      float B = -T(4.0f) * tn * cfg.inv_pow_sigma_z;
      float C = T(2.0f) * cfg.inv_pow_sigma_z;
      float B_2 = (B / T(2.0f)) * (B / T(2.0f));

      float bb_y = bby(A, C, B_2);

      float bb_z = bbz(A, C, B_2);

      int2 center_pixel = pixel_location(y,
                                         z,
                                         cfg.pixel_size,
                                         cfg.pixel_size,
                                         cfg.grid_size_y_,
                                         cfg.grid_size_z_);

      // bounding box limits for event
      int2 ur =
          make_int2(center_pixel.x - pixels_in_line(bb_y, cfg.pixel_size),
                    center_pixel.y + pixels_in_line(bb_z, cfg.pixel_size));
      int2 dl =
          make_int2(center_pixel.x + pixels_in_line(bb_y, cfg.pixel_size),
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
                             tex2D<float>(tex, iy, iz);

            acc += event_kernel * tex2D<float>(tex, iy, iz) *
                   rho[IMAGE_SPACE_LINEAR_INDEX(iy, iz)];
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
                             tex2D<float>(tex, iy, iz);

            atomicAdd(
                &image_buffor[BUFFOR_LINEAR_INDEX(iy, iz)],
                (event_kernel * rho[IMAGE_SPACE_LINEAR_INDEX(iy, iz)]) / acc);
          }
        }
      }
    }
  }
}
