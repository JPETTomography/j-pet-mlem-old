#pragma once

#include <cuda_runtime.h>

#include "../kernel.h"
#include "../strip_detector.h"

#include "config.h"
#include "soa.cuh"
#include "reconstruction_methods.cuh"

template <typename F>
__global__ void reconstruction_2d_strip_cuda(
    StripDetector<F> detector,
    SOA::Events<F>* events,
    int event_list_size,
    F* output_image,
    F* rho,
    cudaTextureObject_t sensitivity_tex,
    int n_blocks,
    int n_threads_per_block) {
  int tid = (blockIdx.x * blockDim.x) + threadIdx.x;

  int block_chunk =
      int(ceilf(event_list_size / (n_blocks * n_threads_per_block)));

  for (int i = 0; i < block_chunk; ++i) {

    if ((i * n_blocks * n_threads_per_block) + tid < event_list_size) {

      for (int j = 0; j < 1; ++j) {

        F y = events->z_u[(i * n_blocks * n_threads_per_block) + tid];
        F z = events->z_d[(i * n_blocks * n_threads_per_block) + tid];
        F acc = 0;

        int y_step = 3 * (detector.sigma_dl / detector.pixel_height);
        int z_step = 3 * (detector.sigma_z / detector.pixel_width);

        Pixel<> center_pixel = detector.pixel_location(y, z);

        for (int iy = center_pixel.x - y_step; iy < center_pixel.x + y_step;
             ++iy) {
          for (int iz = center_pixel.y - z_step; iz < center_pixel.y + z_step;
               ++iz) {

            Point<F> point = detector.pixel_center(iy, iz);
            Kernel<F> kernel;
            float event_kernel =
                kernel.test(y, z, point, detector.sigma_dl, detector.sigma_z);

            acc += event_kernel * rho[IMAGE_SPACE_LINEAR_INDEX(iy, iz)];
          }
        }

        float inv_acc = 1 / acc;

        for (int iz = center_pixel.y - z_step; iz < center_pixel.y + z_step;
             ++iz) {
          for (int iy = center_pixel.x - y_step; iy < center_pixel.x + y_step;
               ++iy) {

            Point<F> point = detector.pixel_center(iy, iz);

            Kernel<F> kernel;
            F event_kernel =
                kernel.test(y, z, point, detector.sigma_dl, detector.sigma_z);

            atomicAdd(&output_image[BUFFER_LINEAR_INDEX(iy, iz)],
                      (event_kernel * rho[IMAGE_SPACE_LINEAR_INDEX(iy, iz)]) *
                          inv_acc);
          }
        }
      }
    }
  }
}
