#pragma once

#include <cuda_runtime.h>

#include "2d/geometry/point.h"
#include "../response.h"
#include "../scanner.h"

#define PIXEL_INDEX(p) (((p).y * scanner.n_z_pixels) + (p).x)

namespace PET2D {
namespace Strip {
namespace GPU {
namespace Reconstruction {

template <template <typename FType> class Kernel, typename F>
__global__ void reconstruction(Scanner<F> scanner,
                               F* responses_z_u,
                               F* responses_z_d,
                               F* responses_dl,
                               const int n_responses,
                               F* output_rho) {
  using Response = Strip::Response<F>;

  Kernel<F> kernel(scanner.sigma_z, scanner.sigma_dl);

  const auto tid = (blockIdx.x * blockDim.x) + threadIdx.x;
  const auto n_threads = gridDim.x * blockDim.x;
  const auto n_chunks = (n_responses + n_threads - 1) / n_threads;

  for (int chunk = 0; chunk < n_chunks; ++chunk) {
    int response_index = chunk * n_threads + tid;
    // check if we are still on the list
    if (response_index >= n_responses) {
      break;
    }

    Response response(responses_z_u[response_index],
                      responses_z_d[response_index],
                      responses_dl[response_index]);

    F tan, y, z;
    response.calculate_tan_y_z(scanner.radius, tan, y, z);

    F sec, A, B, C, bb_y, bb_z;
    kernel.ellipse_bb(tan, sec, A, B, C, bb_y, bb_z);

    Point ellipse_center(z, y);
    Pixel center_pixel = scanner.pixel_at(ellipse_center);

    // bounding box limits for response
    const int bb_half_width = n_pixels_in_line(bb_z, scanner.pixel_width);
    const int bb_half_height = n_pixels_in_line(bb_y, scanner.pixel_height);
    const Pixel bb_tl(center_pixel.x - bb_half_width,
                      center_pixel.y - bb_half_height);
    const Pixel bb_br(center_pixel.x + bb_half_width,
                      center_pixel.y + bb_half_height);

    F denominator = 0;

    for (int iy = bb_tl.y; iy < bb_br.y; ++iy) {
      for (int iz = bb_tl.x; iz < bb_br.x; ++iz) {
        Pixel pixel(iz, iy);
        Point point = scanner.pixel_center(pixel);

        if (kernel.in_ellipse(A, B, C, ellipse_center, point)) {
          point -= ellipse_center;

          F pixel_sensitivity =
              USE_SENSITIVITY ? tex2D(tex_sensitivity, pixel.x, pixel.y) : 1;

          F kernel_value =
              USE_KERNEL ? kernel(y, tan, sec, scanner.radius, point) : 1;

          denominator += kernel_value * tex2D(tex_rho, pixel.x, pixel.y) *
                         pixel_sensitivity;
        }
      }
    }

    F inv_denominator = 1 / denominator;

    for (int iy = bb_tl.y; iy < bb_br.y; ++iy) {
      for (int iz = bb_tl.x; iz < bb_br.x; ++iz) {
        Pixel pixel(iz, iy);
        Point point = scanner.pixel_center(pixel);

        if (kernel.in_ellipse(A, B, C, ellipse_center, point)) {
          point -= ellipse_center;

          F pixel_sensitivity =
              USE_SENSITIVITY ? tex2D(tex_sensitivity, pixel.x, pixel.y) : 1;

          F kernel_value =
              USE_KERNEL ? kernel(y, tan, sec, scanner.radius, point) : 1;

          atomicAdd(&output_rho[PIXEL_INDEX(pixel)],
                    kernel_value * tex2D(tex_rho, pixel.x, pixel.y) /
                        pixel_sensitivity * inv_denominator);
        }
      }
    }
  }
}

template <typename F> _ int n_pixels_in_line(F length, F pixel_size) {
  return (length + F(0.5)) / pixel_size;
}

}  // Reconstruction
}  // GPU
}  // Strip
}  // PET2D
