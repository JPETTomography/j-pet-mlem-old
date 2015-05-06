#pragma once

#include <cuda_runtime.h>

#include "2d/geometry/point.h"
#include "../event.h"
#include "../kernel.h"
#include "../scanner.h"

namespace PET2D {
namespace Strip {
namespace GPU {

template <template <typename Float> class Kernel, typename F>
__global__ void reconstruction(Scanner<F> scanner,
                               F* events_z_u,
                               F* events_z_d,
                               F* events_dl,
                               const int n_events,
                               F* output_rho,
                               const int n_blocks,
                               const int n_threads_per_block) {
  using Point = PET2D::Point<F>;
  using Pixel = PET2D::Pixel<>;
  using Event = Strip::Event<F>;

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

    int y_step = 3 * (scanner.sigma_dl / scanner.pixel_height);
    int z_step = 3 * (scanner.sigma_z / scanner.pixel_width);

    Point ellipse_center(y, z);
    Pixel center_pixel = scanner.pixel_at(ellipse_center);

    for (int iy = center_pixel.y - y_step; iy < center_pixel.y + y_step; ++iy) {
      for (int iz = center_pixel.x - z_step; iz < center_pixel.x + z_step;
           ++iz) {
        Pixel pixel(iz, iy);
        Point point = scanner.pixel_center(pixel);
        Kernel kernel;
        float event_kernel =
            kernel.test(y, z, point, scanner.sigma_dl, scanner.sigma_z);

        denominator += event_kernel * tex2D(tex_rho, pixel.x, pixel.y);
      }
    }

    F inv_denominator = 1 / denominator;

    for (int iy = center_pixel.y - y_step; iy < center_pixel.y + y_step; ++iy) {
      for (int iz = center_pixel.x - z_step; iz < center_pixel.x + z_step;
           ++iz) {
        Pixel pixel(iz, iy);
        Point point = scanner.pixel_center(pixel);

        Kernel kernel;
        F event_kernel =
            kernel.test(y, z, point, scanner.sigma_dl, scanner.sigma_z);

        atomicAdd(
            &output_rho[PIXEL_INDEX(pixel)],
            event_kernel * tex2D(tex_rho, pixel.x, pixel.y) * inv_denominator);
      }
    }
  }
}
}  // GPU
}  // Strip
}  // PET2D
