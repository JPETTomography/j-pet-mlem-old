#pragma once

#include <cuda_runtime.h>

#include "config.h"
#include "event.cuh"
#include "reconstruction_methods.cuh"

template <typename T>
__global__ void reconstruction_2d_strip_cuda(CUDA::Config cfg,
                                             soa_event<float>* soa_data,
                                             Event<T>* event_list,
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

  T sqrt_det_correlation_matrix =
      sqrt(cfg.inv_pow_sigma_z * cfg.inv_pow_sigma_z * cfg.inv_pow_sigma_dl);

  for (int i = 0; i < block_chunk; ++i) {

    if ((i * cfg.number_of_blocks * cfg.number_of_threads_per_block) + tid <
        event_list_size) {

      for (int j = 0; j < 1; ++j) {

#ifdef AOS_ACCESS
        T z_u = event_list[(i * cfg.number_of_blocks *
                            cfg.number_of_threads_per_block) +
                           tid].z_u;
        T z_d = event_list[(i * cfg.number_of_blocks *
                            cfg.number_of_threads_per_block) +
                           tid].z_d;
        T delta_l = event_list[(i * cfg.number_of_blocks *
                                cfg.number_of_threads_per_block) +
                               tid].dl;
#else
        T z_u = soa_data->z_u[(i * cfg.number_of_blocks *
                               cfg.number_of_threads_per_block) +
                              tid];
        T z_d = soa_data->z_d[(i * cfg.number_of_blocks *
                               cfg.number_of_threads_per_block) +
                              tid];
        T delta_l = soa_data->dl[(i * cfg.number_of_blocks *
                                  cfg.number_of_threads_per_block) +
                                 tid];
#endif

        T acc = 0.f;

        T half_grid_size = 0.5f * cfg.grid_size_y;
        T half_pixel_size = 0.5f * cfg.pixel_size;

        // angle space transformation
        T tn = event_tan(z_u, z_d, cfg.R_distance);
        T y = event_y(delta_l, tn);
        T z = event_z(z_u, z_d, y, tn);

        T angle = atanf(tn);

        T cos_ = __cosf(angle);

        T sec_ = float(1.0f) / cos_;
        T sec_sq_ = sec_ * sec_;

        T A = (((T(4.0f) / (cos_ * cos_)) * cfg.inv_pow_sigma_dl) +
               (T(2.0f) * tn * tn * cfg.inv_pow_sigma_z));
        T B = -T(4.0f) * tn * cfg.inv_pow_sigma_z;
        T C = T(2.0f) * cfg.inv_pow_sigma_z;
        T B_2 = (B / T(2.0f)) * (B / T(2.0f));
#if DEBUG
        if (tid == 0 && i == 0) {
          printf("A: %f B: %f C: %f B_2: %f\n", A, B, C, B_2);
        }
#endif
        T bb_y = bby(A, C, B_2);
        T bb_z = bbz(A, C, B_2);
#if DEBUG
        if (tid == 0 && i == 0) {
          printf("bb_y: %f bb_z: %f \n", bb_y, bb_z);
        }
#endif
        int2 center_pixel = pixel_location(y,
                                           z,
                                           cfg.pixel_size,
                                           cfg.pixel_size,
                                           cfg.grid_size_y,
                                           cfg.grid_size_z);

        // bounding box limits for event
        int2 ur =
            make_int2(center_pixel.x - pixels_in_line(bb_y, cfg.pixel_size),
                      center_pixel.y + pixels_in_line(bb_z, cfg.pixel_size));
        int2 dl =
            make_int2(center_pixel.x + pixels_in_line(bb_y, cfg.pixel_size),
                      center_pixel.y - pixels_in_line(bb_z, cfg.pixel_size));
#if DEBUG
        if (tid == 0 && i == 0) {
          printf("UR:= %d Limit:= %d \n", ur.x, ur.y);
          printf("DL:= %d Limit:= %d \n", dl.x, dl.y);
          printf("iz:= %d Limit:= %d \n", dl.y, ur.y);
          printf("iy:= %d Limit:= %d \n", ur.x, dl.x);
        }
#endif
        float2 pp;

        int iter = 0;

        for (int iy = ur.x; iy < dl.x; ++iy) {
          for (int iz = dl.y; iz < ur.y; ++iz) {

            pp = pixel_center(iy,
                              iz,
                              cfg.pixel_size,
                              cfg.pixel_size,
                              cfg.grid_size_y,
                              cfg.grid_size_z,
                              half_grid_size,
                              half_pixel_size);

            if (in_ellipse(A, B, C, y, z, pp)) {
#if DEBUG
              if (tid == 0 && i == 0) {
                printf("Pixel(%d,%d): SUB: %f %f ITER: %d\n",
                       iy,
                       iz,
                       pp.x,
                       pp.y,
                       iter);
              }
#endif
              ++iter;

              pp.x -= y;
              pp.y -= z;

              T event_kernel = main_kernel<T>(y,
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

        if (tid == 0 && i == 0) {

          printf("ACC: %e\n", acc);
        }

        float inv_acc = 1.0f / acc;

        for (int iz = dl.y; iz < ur.y; ++iz) {
          for (int iy = ur.x; iy < dl.x; ++iy) {

            pp = pixel_center(iy,
                              iz,
                              cfg.pixel_size,
                              cfg.pixel_size,
                              cfg.grid_size_y,
                              cfg.grid_size_z,
                              half_grid_size,
                              half_pixel_size);

            if (in_ellipse(A, B, C, y, z, pp)) {

              pp.x -= y;
              pp.y -= z;

              T event_kernel = main_kernel<T>(y,
                                              tn,
                                              sec_,
                                              sec_sq_,
                                              pp,
                                              inv_c,
                                              cfg,
                                              sqrt_det_correlation_matrix) /
                               tex2D<float>(tex, iy, iz);

              atomicAdd(&image_buffor[BUFFOR_LINEAR_INDEX(iy, iz)],
                        (event_kernel * rho[IMAGE_SPACE_LINEAR_INDEX(iy, iz)]) *
                            inv_acc);
            }
          }
        }
      }
    }
  }
}
