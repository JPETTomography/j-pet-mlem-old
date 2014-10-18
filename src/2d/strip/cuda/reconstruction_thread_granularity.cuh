#pragma once

#include <cuda_runtime.h>

#include "2d/geometry/point.h"
#include "../event.h"
#include "../detector.h"

#define PIXEL_INDEX(p) (((p).y * detector.n_z_pixels) + (p).x)

namespace PET2D {
namespace Strip {
namespace GPU {

template <template <typename Float> class Kernel, typename F>
__global__ void reconstruction(Detector<F> detector,
                               F* events_z_u,
                               F* events_z_d,
                               F* events_dl,
                               const int n_events,
                               F* output_rho,
                               const int n_blocks,
                               const int n_threads_per_block) {
  Kernel<F> kernel(detector.sigma_z, detector.sigma_dl);

  int n_threads = n_blocks * n_threads_per_block;
  int n_chunks = (n_events + n_threads - 1) / n_threads;

  int tid = (blockIdx.x * blockDim.x) + threadIdx.x;

  for (int chunk = 0; chunk < n_chunks; ++chunk) {
    int event_index = chunk * n_threads + tid;
    // check if we are still on the list
    if (event_index >= n_events) {
      break;
    }

    Event<F> event(events_z_u[event_index],
                   events_z_d[event_index],
                   events_dl[event_index]);

    F tan, y, z;
    event.transform(detector.radius, tan, y, z);

    F sec, A, B, C, bb_y, bb_z;
    kernel.ellipse_bb(tan, sec, A, B, C, bb_y, bb_z);

    Point<F> ellipse_center(z, y);
    Pixel<> center_pixel = detector.pixel_location(ellipse_center);

    // bounding box limits for event
    const int bb_half_width = n_pixels_in_line(bb_z, detector.pixel_width);
    const int bb_half_height = n_pixels_in_line(bb_y, detector.pixel_height);
    const Pixel<> tl(center_pixel.x - bb_half_width,
                     center_pixel.y - bb_half_height);
    const Pixel<> br(center_pixel.x + bb_half_width,
                     center_pixel.y + bb_half_height);

    F denominator = 0;

    for (int iy = tl.y; iy < br.y; ++iy) {
      for (int iz = tl.x; iz < br.x; ++iz) {
        Pixel<> pixel(iz, iy);
        Point<F> point = detector.pixel_center(pixel);

        if (kernel.in_ellipse(A, B, C, ellipse_center, point)) {
          point -= ellipse_center;

          F pixel_sensitivity =
              USE_SENSITIVITY ? tex2D(tex_sensitivity, pixel.x, pixel.y) : 1;

          F event_kernel =
              USE_KERNEL ? kernel(y, tan, sec, detector.radius, point) : 1;

          denominator += event_kernel * tex2D(tex_rho, pixel.x, pixel.y) *
                         pixel_sensitivity;
        }
      }
    }

    F inv_denominator = 1 / denominator;

    for (int iy = tl.y; iy < br.y; ++iy) {
      for (int iz = tl.x; iz < br.x; ++iz) {
        Pixel<> pixel(iz, iy);
        Point<F> point = detector.pixel_center(pixel);

        if (kernel.in_ellipse(A, B, C, ellipse_center, point)) {
          point -= ellipse_center;

          F pixel_sensitivity =
              USE_SENSITIVITY ? tex2D(tex_sensitivity, pixel.x, pixel.y) : 1;

          F event_kernel =
              USE_KERNEL ? kernel(y, tan, sec, detector.radius, point) : 1;

          atomicAdd(&output_rho[PIXEL_INDEX(pixel)],
                    event_kernel * tex2D(tex_rho, pixel.x, pixel.y) /
                        pixel_sensitivity * inv_denominator);
        }
      }
    }
  }
}

template <typename F> _ int n_pixels_in_line(F length, F pixel_size) {
  return (length + F(0.5)) / pixel_size;
}
}  // GPU
}  // Strip
}  // PET2D
