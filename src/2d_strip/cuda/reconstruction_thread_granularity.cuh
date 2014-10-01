#pragma once

#include <cuda_runtime.h>

#include "geometry/point.h"
#include "../event.h"
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
  Kernel<F> kernel;

#if USE_RUNTIME_TRUE_FALSE
  // The variables below always evaluate to true/false, but compiler cannot
  // assume that since n_blocks is only known at runtime.
  bool rt_true = /***/ (n_blocks > 0);
  bool rt_false = /**/ (n_blocks == 0);
#endif

  // In optimized build we set it to true/false which triggers else branches to
  // be optimized out of the code, however in code parts benchmark we shall use
  // rt_true/rt_false that guarantees else branches to be not optimized, so we
  // can reliably measure time disabling certain computations.
  bool use_kernel = true;
  bool use_sensitivity = false;

  F sqrt_det_cor_mat = detector.sqrt_det_cor_mat();
  int n_threads = n_blocks * n_threads_per_block;
  int n_chunks = (n_events + n_threads - 1) / n_threads;

  int tid = (blockIdx.x * blockDim.x) + threadIdx.x;

  for (int chunk = 0; chunk < n_chunks; ++chunk) {
    int i = chunk * n_threads + tid;
    // check if we are still on the list
    if (i >= n_events) {
      break;
    }

#if AOS_ACCESS
    Event<F> event = events[i];
#else
    Event<F> event(events_z_u[i], events_z_d[i], events_dl[i]);
#endif

    F tan, y, z;
    event.transform(detector.radius, tan, y, z);
    F angle = compat::atan(tan);
    Point<F> ellipse_center(z, y);

    F sec, sec_sq, A, B, C, bb_y, bb_z;
    detector.ellipse_bb(angle, tan, sec, sec_sq, A, B, C, bb_y, bb_z);

    Pixel<> center_pixel = detector.pixel_location(ellipse_center);

    // bounding box limits for event
    const int bb_half_width = n_pixels_in_line(bb_z, detector.pixel_width);
    const int bb_half_height = n_pixels_in_line(bb_y, detector.pixel_height);
    const Pixel<> tl(center_pixel.x - bb_half_width,
                     center_pixel.y - bb_half_height);
    const Pixel<> br(center_pixel.x + bb_half_width,
                     center_pixel.y + bb_half_height);

    F acc = 0;

    for (int iy = tl.y; iy < br.y; ++iy) {
      for (int iz = tl.x; iz < br.x; ++iz) {
        Pixel<> pixel(iz, iy);
        Point<F> point = detector.pixel_center(pixel);

        if (detector.in_ellipse(A, B, C, ellipse_center, point)) {
          point -= ellipse_center;

#if USE_KERNEL
          F event_kernel = kernel(y,
                                  tan,
                                  sec,
                                  sec_sq,
                                  detector.radius,
                                  point,
                                  detector.inv_cor_mat_diag,
                                  sqrt_det_cor_mat) *
                           inv_pixel_sensitivity;
#else
          F event_kernel = 1;
#endif

          acc += event_kernel * tex2D(tex_rho, pixel.x, pixel.y);
        }
      }
    }

    F inv_acc = 1 / acc;

    for (int iy = tl.y; iy < br.y; ++iy) {
      for (int iz = tl.x; iz < br.x; ++iz) {
        Pixel<> pixel(iz, iy);
        Point<F> point = detector.pixel_center(pixel);

        if (detector.in_ellipse(A, B, C, ellipse_center, point)) {
          point -= ellipse_center;

          F pixel_sensitivity =
              use_sensitivity ? tex2D(tex_sensitivity, pixel.x, pixel.y) : 1;

          F event_kernel = use_kernel ? kernel(y,
                                               tan,
                                               sec,
                                               sec_sq,
                                               detector.radius,
                                               point,
                                               detector.inv_cor_mat_diag,
                                               sqrt_det_cor_mat)
                                      : 1;

          atomicAdd(&output_rho[PIXEL_INDEX(pixel)],
                    event_kernel * tex2D(tex_rho, pixel.x, pixel.y) /
                        pixel_sensitivity * inv_acc);
        }
      }
    }
  }
}

template <typename F> _ int n_pixels_in_line(F length, F pixel_size) {
  return (length + F(0.5)) / pixel_size;
}
