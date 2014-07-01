#pragma once

#include <cuda_runtime.h>

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

  int block_chunk =
      int(ceilf(event_list_size / (cfg.n_blocks * cfg.n_threads_per_block)));

  for (int i = 0; i < block_chunk; ++i) {

    if ((i * cfg.n_blocks * cfg.n_threads_per_block) + tid < event_list_size) {

      for (int j = 0; j < 1; ++j) {

        F y = events->z_u[(i * cfg.n_blocks * cfg.n_threads_per_block) + tid];
        F z = events->z_d[(i * cfg.n_blocks * cfg.n_threads_per_block) + tid];
        F acc = 0;

        if (tid == 0 && i == 0) {

          printf("%f %f\n", y, z);
        }

        F half_grid_size = 0.5f * cfg.grid_size_y;
        F half_pixel_size = 0.5f * cfg.pixel_size;

        int y_step = 3 * (cfg.dl / cfg.pixel_size);
        int z_step = 3 * (cfg.sigma / cfg.pixel_size);

        Pixel<> center_pixel = pixel_location(y,
                                              z,
                                              cfg.pixel_size,
                                              cfg.pixel_size,
                                              cfg.grid_size_y,
                                              cfg.grid_size_z);
#if DEBUG
        if (tid == 0 && i == 0) {
          printf("TID: %d %f %f LIMIT: %d%d\n", tid, y, z, y_step, z_step);
          printf("TID: %d %d %d %d %d\n",
                 tid,
                 center_pixel.x - z_step,
                 center_pixel.x + z_step,
                 center_pixel.y - y_step,
                 center_pixel.y + y_step);
        }
#endif
        for (int iy = center_pixel.x - y_step; iy < center_pixel.x + y_step;
             ++iy) {
          for (int iz = center_pixel.y - z_step; iz < center_pixel.y + z_step;
               ++iz) {

            Point<float> point = pixel_center(iy,
                                              iz,
                                              cfg.pixel_size,
                                              cfg.pixel_size,
                                              half_grid_size,
                                              half_pixel_size);

            float event_kernel = test_kernel<float>(y, z, point, cfg);

            acc += event_kernel * rho[IMAGE_SPACE_LINEAR_INDEX(iy, iz)];

#if DEBUG
            if (tid == 0 && i == 0) {
              printf("PP: %f %f %e EVENT:  ", point.x, point.y, event_kernel);
              printf("TID: %d %d %d \n", tid, iy, iz);
              printf("INDEX: %d\n", BUFFOR_LINEAR_INDEX(iy, iz));
            }
#endif
          }
        }

        float inv_acc = 1 / acc;

        for (int iz = center_pixel.y - z_step; iz < center_pixel.y + z_step;
             ++iz) {
          for (int iy = center_pixel.x - y_step; iy < center_pixel.x + y_step;
               ++iy) {

            Point<F> point = pixel_center(iy,
                                          iz,
                                          cfg.pixel_size,
                                          cfg.pixel_size,
                                          half_grid_size,
                                          half_pixel_size);

            F event_kernel = test_kernel<F>(y, z, point, cfg);

            atomicAdd(&output_image[BUFFER_LINEAR_INDEX(iy, iz)],
                      (event_kernel * rho[IMAGE_SPACE_LINEAR_INDEX(iy, iz)]) *
                          inv_acc);
          }
        }
      }
    }
  }
}
