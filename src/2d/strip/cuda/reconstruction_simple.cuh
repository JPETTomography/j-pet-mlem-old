#pragma once

#include <cuda_runtime.h>

#include "2d/geometry/point.h"
#include "../response.h"
#include "../kernel.h"
#include "../scanner.h"

namespace PET2D {
namespace Strip {
namespace GPU {

template <template <typename Float> class Kernel, typename F>
__global__ void reconstruction(Scanner<F> scanner,
                               F* responses_z_u,
                               F* responses_z_d,
                               F* responses_dl,
                               const int n_responses,
                               F* output_rho) {
  using Point = PET2D::Point<F>;
  using Pixel = PET2D::Pixel<>;
  using Response = Strip::Response<F>;

  // mark variables used
  (void)(responses_dl);

  const auto tid = (blockIdx.x * blockDim.x) + threadIdx.x;
  const auto n_threads = gridDim.x * blockDim.x;
  const auto n_chunks = (n_responses + n_threads - 1) / n_threads;

  for (int chunk = 0; chunk < n_chunks; ++chunk) {
    int i = chunk * n_threads + tid;
    // check if we are still on the list
    if (i >= n_responses) {
      break;
    }

    F y, z;
    y = responses_z_u[i];
    z = responses_z_d[i];
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
        F kernel_value =
            kernel.test(y, z, point, scanner.sigma_dl, scanner.sigma_z);

        denominator += kernel_value * tex2D(tex_rho, pixel.x, pixel.y);
      }
    }

    F inv_denominator = 1 / denominator;

    for (int iy = center_pixel.y - y_step; iy < center_pixel.y + y_step; ++iy) {
      for (int iz = center_pixel.x - z_step; iz < center_pixel.x + z_step;
           ++iz) {
        Pixel pixel(iz, iy);
        Point point = scanner.pixel_center(pixel);

        Kernel kernel;
        F kernel_value =
            kernel.test(y, z, point, scanner.sigma_dl, scanner.sigma_z);

        atomicAdd(
            &output_rho[PIXEL_INDEX(pixel)],
            kernel_value * tex2D(tex_rho, pixel.x, pixel.y) * inv_denominator);
      }
    }
  }
}

}  // GPU
}  // Strip
}  // PET2D
