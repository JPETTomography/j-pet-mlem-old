#pragma once

#include <cuda_runtime.h>

#include "geometry/point.h"
#include "../event.h"
#include "../kernel.h"
#include "../strip_detector.h"

#include "config.h"

template <typename F> _ int n_pixels_in_line(F length, F pixel_size) {
  return (length + F(0.5)) / pixel_size;
}

template <typename F>
__global__ void reconstruction(StripDetector<F> detector,
                               F* events_z_u,
                               F* events_z_d,
                               F* events_dl,
                               int n_events,
                               F* output,
                               F* rho,
                               TEX_ARG(sensitivity),
                               int n_blocks,
                               int n_threads_per_block) {
  Kernel<F> kernel;
  F sqrt_det_cor_mat = detector.sqrt_det_cor_mat();

  int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
  int n_threads = n_blocks * n_threads_per_block;
  int n_chunks = (n_events + n_threads - 1) / n_threads;

  for (int chunk = 0; chunk < n_chunks; ++chunk) {
    int i = chunk * n_threads + tid;
    // check if we are still on the list
    if (i >= n_events) {
      break;
    }

#ifdef AOS_ACCESS
    Event<F> event = events[i];
#else
    Event<F> event(events_z_u[i], events_z_d[i], events_dl[i]);
#endif

    F tan, y, z;
    event.transform(detector.radius, tan, y, z);
    F angle = compat::atan(tan);
    Point<F> ellipse_center(y, z);

    F sec, sec_sq, A, B, C, bb_y, bb_z;
    detector.ellipse_bb(angle, tan, sec, sec_sq, A, B, C, bb_y, bb_z);

    Pixel<> center_pixel = detector.pixel_location(y, z);

    // bounding box limits for event
    Pixel<> ur(center_pixel.x - n_pixels_in_line(bb_y, detector.pixel_height),
               center_pixel.y + n_pixels_in_line(bb_z, detector.pixel_width));
    Pixel<> dl(center_pixel.x + n_pixels_in_line(bb_y, detector.pixel_height),
               center_pixel.y - n_pixels_in_line(bb_z, detector.pixel_width));

    F acc = 0;

    for (int iy = ur.x; iy < dl.x; ++iy) {
      for (int iz = dl.y; iz < ur.y; ++iz) {
        Point<F> point = detector.pixel_center(iy, iz);

        if (detector.in_ellipse(A, B, C, ellipse_center, point)) {
          point -= ellipse_center;

          F event_kernel = kernel(y,
                                  tan,
                                  sec,
                                  sec_sq,
                                  detector.radius,
                                  point,
                                  detector.inv_cor_mat_diag,
                                  sqrt_det_cor_mat) /
                           TEX_2D(F, sensitivity, iy, iz);

          acc += event_kernel * TEX_2D(F, sensitivity, iy, iz) *
                 rho[IMAGE_SPACE_LINEAR_INDEX(iy, iz)];
        }
      }
    }

    F inv_acc = 1 / acc;

    for (int iz = dl.y; iz < ur.y; ++iz) {
      for (int iy = ur.x; iy < dl.x; ++iy) {
        Point<F> point = detector.pixel_center(iy, iz);

        if (detector.in_ellipse(A, B, C, ellipse_center, point)) {
          point -= ellipse_center;

          F event_kernel = kernel(y,
                                  tan,
                                  sec,
                                  sec_sq,
                                  detector.radius,
                                  point,
                                  detector.inv_cor_mat_diag,
                                  sqrt_det_cor_mat) /
                           TEX_2D(F, sensitivity, iy, iz);

          atomicAdd(
              &output[BUFFER_LINEAR_INDEX(iy, iz)],
              event_kernel * rho[IMAGE_SPACE_LINEAR_INDEX(iy, iz)] * inv_acc);
        }
      }
    }
  }
}
