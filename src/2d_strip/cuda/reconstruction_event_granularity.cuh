#pragma once

#include <cuda_runtime.h>

#include "geometry/point.h"
#include "../event.h"
#include "../kernel.h"
#include "../strip_detector.h"

#include "config.h"
#include "soa.cuh"
#include "reconstruction_methods.cuh"

template <typename F>
__global__ void reconstruction_2d_strip_cuda(StripDetector<F> detector,
                                             SOA::Events<F>* events,
                                             int n_events,
                                             F* output_image,
                                             F* rho,
                                             cudaTextureObject_t sensitivity,
                                             int n_blocks,
                                             int n_threads_per_block) {

  Kernel<F> kernel;
  F sqrt_det_cor_mat = detector.sqrt_det_cor_mat();

  int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
  int block_chunk = int(ceilf(n_events / (n_blocks * n_threads_per_block)));

  for (int i = 0; i < block_chunk; ++i) {
    // check if we are still on the list
    if (i * n_blocks * n_threads_per_block + tid >= n_events) {
      break;
    }

    for (int j = 0; j < 1; ++j) {

#ifdef AOS_ACCESS
      Event<F> event =
          events[(i * cfg.number_of_blocks * cfg.number_of_threads_per_block) +
                 tid];
#else
      Event<F> event(events->z_u[(i * n_blocks * n_threads_per_block) + tid],
                     events->z_d[(i * n_blocks * n_threads_per_block) + tid],
                     events->dl[(i * n_blocks * n_threads_per_block) + tid]);
#endif

      F acc = 0;

      // angle space transformation
      F tan = event.tan(detector.radius);
      F y = event.y(tan);
      F z = event.z(y, tan);
      F angle = compat::atan(tan);
      Point<F> ellipse_center(y, z);

      F sec, sec_sq, A, B, C, bb_y, bb_z;
      detector.ellipse_bb(angle, tan, sec, sec_sq, A, B, C, bb_y, bb_z);

      Pixel<> center_pixel = detector.pixel_location(y, z);

      // bounding box limits for event
      Pixel<> ur(center_pixel.x - pixels_in_line(bb_y, detector.pixel_height),
                 center_pixel.y + pixels_in_line(bb_z, detector.pixel_width));
      Pixel<> dl(center_pixel.x + pixels_in_line(bb_y, detector.pixel_height),
                 center_pixel.y - pixels_in_line(bb_z, detector.pixel_width));

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
                             tex2D<F>(sensitivity, iy, iz);

            acc += event_kernel * tex2D<F>(sensitivity, iy, iz) *
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
                             tex2D<F>(sensitivity, iy, iz);

            atomicAdd(&output_image[BUFFER_LINEAR_INDEX(iy, iz)],
                      (event_kernel * rho[IMAGE_SPACE_LINEAR_INDEX(iy, iz)]) *
                          inv_acc);
          }
        }
      }
    }
  }
}
