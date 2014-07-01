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
    cudaTextureObject_t sensitiviy_tex) {
  short sh_mem_pixel_buffor[20 * 512];

  int tid = (blockIdx.x * blockDim.x) + threadIdx.x;

  F inv_c[3] = { cfg.inv_pow_sigma_z,
                 cfg.inv_pow_sigma_z,
                 cfg.inv_pow_sigma_dl };

  F half_grid_size = F(0.5) * cfg.grid_size_y;
  F half_pixel_size = F(0.5) * cfg.pixel_size;

  int offset;

  int number_of_blocks =
      int(ceilf(event_list_size /
                (cfg.n_blocks * (cfg.n_threads_per_block / WARP_SIZE))));

  int block_size = (cfg.n_blocks * (cfg.n_threads_per_block / WARP_SIZE));

  int thread_warp_index = floorf(threadIdx.x / WARP_SIZE);

  F sqrt_det_correlation_matrix =
      sqrt(cfg.inv_pow_sigma_z * cfg.inv_pow_sigma_z * cfg.inv_pow_sigma_dl);

  for (int i = 0; i < number_of_blocks; ++i) {

    // calculate offset in event memory location for warp
    int warp_id =
        (i * block_size + (blockIdx.x * (cfg.n_threads_per_block / WARP_SIZE)) +
         thread_warp_index);

    if ((warp_id < event_list_size)) {

      Event<F> event(
          events->z_u[warp_id], events->z_d[warp_id], events->dl[warp_id]);

      F acc = 0;

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
      Pixel<> ul(center_pixel.x - pixels_in_line(bb_y, cfg.pixel_size),
                 center_pixel.y - pixels_in_line(bb_z, cfg.pixel_size));

#if DEBUG
      if (tid == 0 && i == 0) {
        printf("UR:= %d Limit:= %d \n", ur.x, ur.y);
        printf("DL:= %d Limit:= %d \n", dl.x, dl.y);
        printf("iz:= %d Limit:= %d \n", dl.y, ur.y);
        printf("iy:= %d Limit:= %d \n", ur.x, dl.x);
      }
#endif

      int bb_size = ceilf(((ur.y - ul.y) * (dl.x - ur.x)) / WARP_SIZE);

      int pixel_count = 0;

      Pixel<> first_pixel = ul;

      for (int k = 0; k < bb_size; ++k) {
        Pixel<> pixel = warp_space_pixel(offset, first_pixel, ul, ur, dl, tid);
        Point<F> point = pixel_center(pixel.x,
                                      pixel.y,
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
                           tex2D<F>(sensitiviy_tex, pixel.x, pixel.y);
#if DEBUG
          if (tid < 32 && i == 0) {
                        printf("TID: %d KERNEL: %e PIXEL: %d %d
                        K:%d\n",tid,event_kernel,tid_pixel.x,tid_pixel.y,k);
          }
#endif
          acc += event_kernel * tex2D<F>(sensitiviy_tex, pixel.x, pixel.y) *
                 rho[IMAGE_SPACE_LINEAR_INDEX(pixel.x, pixel.y)];

          sh_mem_pixel_buffor[SH_MEM_INDEX(threadIdx.x, pixel_count, 0)] =
              pixel.x;
          sh_mem_pixel_buffor[SH_MEM_INDEX(threadIdx.x, pixel_count, 1)] =
              pixel.y;
          pixel_count++;
        }
      }

      for (int xor_iter = 16; xor_iter >= 1; xor_iter /= 2) {
        acc += __shfl_xor(acc, xor_iter, WARP_SIZE);
      }

#if DEBUG
      if (tid == 0 && i == 0) {
        printf("KERNEL: %e\n", acc);
      }
#endif

      F inv_acc = 1 / acc;

      for (int k = 0; k < pixel_count; ++k) {

        Pixel<> pixel(sh_mem_pixel_buffor[SH_MEM_INDEX(threadIdx.x, k, 0)],
                      sh_mem_pixel_buffor[SH_MEM_INDEX(threadIdx.x, k, 1)]);
        Point<F> point = pixel_center(pixel.x,
                                      pixel.y,
                                      cfg.pixel_size,
                                      cfg.pixel_size,
                                      half_grid_size,
                                      half_pixel_size);
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
                         tex2D<F>(sensitiviy_tex, pixel.x, pixel.y);

#if DEBUG
        if (tid < 32 && i == 0 && k < loop_index) {
                    printf("TID: %d KERNEL: %e PIXEL: %d %d
                    K:%d\n",tid,event_kernel,tid_pixel.x,tid_pixel.y,k);
        }
#endif

        atomicAdd(
            &output_image[BUFFER_LINEAR_INDEX(pixel.x, pixel.y)],
            (event_kernel * rho[IMAGE_SPACE_LINEAR_INDEX(pixel.x, pixel.y)]) *
                inv_acc * sec_sq);
      }
    }
  }
}
