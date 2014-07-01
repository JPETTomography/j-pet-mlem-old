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

  int block_chunk =
      int(ceilf(event_list_size / (cfg.n_blocks * cfg.n_threads_per_block)));

  for (int i = 0; i < block_chunk; ++i) {

    if ((i * cfg.n_blocks * cfg.n_threads_per_block) + tid < event_list_size) {

      for (int j = 0; j < 1; ++j) {

        T y = soa_data->z_u[(i * cfg.n_blocks * cfg.n_threads_per_block) + tid];
        T z = soa_data->z_d[(i * cfg.n_blocks * cfg.n_threads_per_block) + tid];
        T acc = 0;

        if (tid == 0 && i == 0) {

          printf("%f %f\n", y, z);
        }

        T half_grid_size = 0.5f * cfg.grid_size_y;
        T half_pixel_size = 0.5f * cfg.pixel_size;

        int y_step = 3 * (cfg.dl / cfg.pixel_size);
        int z_step = 3 * (cfg.sigma / cfg.pixel_size);

        int2 center_pixel = pixel_location(y,
                                           z,
                                           cfg.pixel_size,
                                           cfg.pixel_size,
                                           cfg.grid_size_y,
                                           cfg.grid_size_z);
        float2 pp;
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

            pp = pixel_center(iy,
                              iz,
                              cfg.pixel_size,
                              cfg.pixel_size,
                              cfg.grid_size_y,
                              cfg.grid_size_z,
                              half_grid_size,
                              half_pixel_size);

            T event_kernel = test_kernel<T>(y, z, pp, cfg);

            acc += event_kernel * rho[IMAGE_SPACE_LINEAR_INDEX(iy, iz)];

#if DEBUG
            if (tid == 0 && i == 0) {
              printf("PP: %f %f %e EVENT:  ", pp.x, pp.y, event_kernel);
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

            pp = pixel_center(iy,
                              iz,
                              cfg.pixel_size,
                              cfg.pixel_size,
                              cfg.grid_size_y,
                              cfg.grid_size_z,
                              half_grid_size,
                              half_pixel_size);

            T event_kernel = test_kernel<T>(y, z, pp, cfg);

            atomicAdd(&image_buffor[BUFFOR_LINEAR_INDEX(iy, iz)],
                      (event_kernel * rho[IMAGE_SPACE_LINEAR_INDEX(iy, iz)]) *
                          inv_acc);
          }
        }
      }
    }
  }
}
