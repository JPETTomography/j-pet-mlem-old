#pragma once

#include <cuda_runtime.h>

#include "geometry/point.h"
#include "../event.h"
#include "../kernel.h"

#include "config.h"
#include "soa.cuh"
#include "reconstruction_methods.cuh"

template <typename F>
__global__ void reconstruction_2d_strip_cuda(
    CUDA::Config<F> cfg,
    SOA::Events<F>* events,
    int event_list_size,
    F* output_image,
    F* rho,
    cudaTextureObject_t sensitivity_tex) {
  int tid = (blockIdx.x * blockDim.x) + threadIdx.x;

  __shared__ F inv_c[3];

  int block_chunk =
      int(ceilf(event_list_size / (cfg.n_blocks * cfg.n_threads_per_block)));

  if (threadIdx.x == 0) {
    inv_c[0] = cfg.inv_pow_sigma_z;
    inv_c[1] = cfg.inv_pow_sigma_z;
    inv_c[2] = cfg.inv_pow_sigma_dl;
  }

  __syncthreads();

  F sqrt_det_correlation_matrix =
      sqrt(cfg.inv_pow_sigma_z * cfg.inv_pow_sigma_z * cfg.inv_pow_sigma_dl);

  for (int i = 0; i < block_chunk; ++i) {
    // check if we are still on the list
    if (i * cfg.n_blocks * cfg.n_threads_per_block + tid >= event_list_size) {
      break;
    }

    for (int j = 0; j < 1; ++j) {

#ifdef AOS_ACCESS
      T z_u =
          events[(i * cfg.number_of_blocks * cfg.number_of_threads_per_block) +
                 tid].z_u;
      T z_d =
          events[(i * cfg.number_of_blocks * cfg.number_of_threads_per_block) +
                 tid].z_d;
      T delta_l =
          events[(i * cfg.number_of_blocks * cfg.number_of_threads_per_block) +
                 tid].dl;
#else
      Event<F> event(
          events->z_u[(i * cfg.n_blocks * cfg.n_threads_per_block) + tid],
          events->z_d[(i * cfg.n_blocks * cfg.n_threads_per_block) + tid],
          events->dl[(i * cfg.n_blocks * cfg.n_threads_per_block) + tid]);
#endif

      F acc = 0;

      F half_grid_size = F(0.5) * cfg.grid_size_y;
      F half_pixel_size = F(0.5) * cfg.pixel_size;

      // angle space transformation
      F tan = event.tan(cfg.R_distance);
      F y = event.y(tan);
      F z = event.z(y, tan);
      F angle = compat::atan(tan);
      F cos = compat::cos(angle);

      F sec = 1 / cos;
      F sec_sq = sec * sec;

      F A = (((4 / (cos * cos)) * cfg.inv_pow_sigma_dl) +
             (2 * tan * tan * cfg.inv_pow_sigma_z));
      F B = -4 * tan * cfg.inv_pow_sigma_z;
      F C = 2 * cfg.inv_pow_sigma_z;
      F B_2 = (B / 2) * (B / 2);

      F bb_y = bby(A, C, B_2);
      F bb_z = bbz(A, C, B_2);

#if DEBUG
      if (tid == 0 && i == 0) {
        printf("A: %f B: %f C: %f B_2: %f\n", A, B, C, B_2);
        printf("bb_y: %f bb_z: %f \n", bb_y, bb_z);
      }
#endif

      Pixel<> center_pixel = pixel_location(y,
                                            z,
                                            cfg.pixel_size,
                                            cfg.pixel_size,
                                            cfg.grid_size_y,
                                            cfg.grid_size_z);

      // bounding box limits for event
      Pixel<> ur(center_pixel.x - pixels_in_line(bb_y, cfg.pixel_size),
                 center_pixel.y + pixels_in_line(bb_z, cfg.pixel_size));
      Pixel<> dl(center_pixel.x + pixels_in_line(bb_y, cfg.pixel_size),
                 center_pixel.y - pixels_in_line(bb_z, cfg.pixel_size));
#if DEBUG
      if (tid == 0 && i == 0) {
        printf("UR:= %d Limit:= %d \n", ur.x, ur.y);
        printf("DL:= %d Limit:= %d \n", dl.x, dl.y);
        printf("iz:= %d Limit:= %d \n", dl.y, ur.y);
        printf("iy:= %d Limit:= %d \n", ur.x, dl.x);
      }
#endif
      for (int iy = ur.x; iy < dl.x; ++iy) {
        for (int iz = dl.y; iz < ur.y; ++iz) {

          Point<F> point = pixel_center(iy,
                                        iz,
                                        cfg.pixel_size,
                                        cfg.pixel_size,
                                        half_grid_size,
                                        half_pixel_size);

          if (in_ellipse(A, B, C, y, z, point)) {
            point.x -= y;
            point.y -= z;

            Kernel<F> kernel;
            F event_kernel = kernel(y,
                                    tan,
                                    sec,
                                    sec_sq,
                                    cfg.R_distance,
                                    point,
                                    inv_c,
                                    sqrt_det_correlation_matrix) /
                             tex2D<F>(sensitivity_tex, iy, iz);

            acc += event_kernel * tex2D<F>(sensitivity_tex, iy, iz) *
                   rho[IMAGE_SPACE_LINEAR_INDEX(iy, iz)];
          }
        }
      }

      if (tid == 0 && i == 0) {
        printf("ACC: %e\n", acc);
      }

      F inv_acc = 1 / acc;

      for (int iz = dl.y; iz < ur.y; ++iz) {
        for (int iy = ur.x; iy < dl.x; ++iy) {

          Point<F> point = pixel_center(iy,
                                        iz,
                                        cfg.pixel_size,
                                        cfg.pixel_size,
                                        half_grid_size,
                                        half_pixel_size);

          if (in_ellipse(A, B, C, y, z, point)) {

            point.x -= y;
            point.y -= z;

            Kernel<F> kernel;
            F event_kernel = kernel(y,
                                    tan,
                                    sec,
                                    sec_sq,
                                    cfg.R_distance,
                                    point,
                                    inv_c,
                                    sqrt_det_correlation_matrix) /
                             tex2D<F>(sensitivity_tex, iy, iz);

            atomicAdd(&output_image[BUFFER_LINEAR_INDEX(iy, iz)],
                      (event_kernel * rho[IMAGE_SPACE_LINEAR_INDEX(iy, iz)]) *
                          inv_acc);
          }
        }
      }
    }
  }
}
