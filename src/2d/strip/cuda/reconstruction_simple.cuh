#pragma once

#include <cuda_runtime.h>

#include "geometry/point.h"
#include "../event.h"
#include "../kernel.h"
#include "../strip_detector.h"

#include "config.h"

template <template <typename Float> class Kernel, typename F>
__global__ void reconstruction(StripDetector<F> detector,
                               F* events_z_u,
                               F* events_z_d,
                               F* events_dl,
                               const int n_events,
                               F* output_rho,
                               const int n_blocks,
                               const int n_threads_per_block) {
  // mark variables used
  (void)(events_dl);

  int n_threads = n_blocks * n_threads_per_block;
  int n_chunks = (n_events + n_threads - 1) / n_threads;

  int tid = (blockIdx.x * blockDim.x) + threadIdx.x;

  for (int chunk = 0; chunk < n_chunks; ++chunk) {
    int i = chunk * n_threads + tid;
    // check if we are still on the list
    if (i >= n_events) {
      break;
    }

    F y, z;
    y = events_z_u[i];
    z = events_z_d[i];
    F denominator = 0;

    int y_step = 3 * (detector.sigma_dl / detector.pixel_height);
    int z_step = 3 * (detector.sigma_z / detector.pixel_width);

    Point<F> ellipse_center(y, z);
    Pixel<> center_pixel = detector.pixel_location(ellipse_center);

    for (int iy = center_pixel.y - y_step; iy < center_pixel.y + y_step; ++iy) {
      for (int iz = center_pixel.x - z_step; iz < center_pixel.x + z_step;
           ++iz) {
        Pixel<> pixel(iz, iy);
        Point<F> point = detector.pixel_center(pixel);
        Kernel<F> kernel;
        float event_kernel =
            kernel.test(y, z, point, detector.sigma_dl, detector.sigma_z);

        denominator += event_kernel * tex2D(tex_rho, pixel.x, pixel.y);
      }
    }

    F inv_denominator = 1 / denominator;

    for (int iy = center_pixel.y - y_step; iy < center_pixel.y + y_step; ++iy) {
      for (int iz = center_pixel.x - z_step; iz < center_pixel.x + z_step;
           ++iz) {
        Pixel<> pixel(iz, iy);
        Point<F> point = detector.pixel_center(pixel);

        Kernel<F> kernel;
        F event_kernel =
            kernel.test(y, z, point, detector.sigma_dl, detector.sigma_z);

        atomicAdd(
            &output_rho[PIXEL_INDEX(pixel)],
            event_kernel * tex2D(tex_rho, pixel.x, pixel.y) * inv_denominator);
      }
    }
  }
}
