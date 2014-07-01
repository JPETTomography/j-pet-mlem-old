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

  short sh_mem_pixel_buffer[20 * 512];

  int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
  int offset;
  int block_size = (n_blocks * (n_threads_per_block / WARP_SIZE));
  int number_of_blocks = (n_events + block_size - 1) / block_size;
  int thread_warp_index = threadIdx.x / WARP_SIZE;

  for (int i = 0; i < number_of_blocks; ++i) {

    // calculate offset in event memory location for warp
    int warp_id =
        (i * block_size + (blockIdx.x * (n_threads_per_block / WARP_SIZE)) +
         thread_warp_index);

    if (warp_id >= n_events)
      break;

    Event<F> event(
        events->z_u[warp_id], events->z_d[warp_id], events->dl[warp_id]);

    F acc = 0;

    // angle space transformation
    F tan = event.tan(detector.radius);
    F y = event.y(tan);
    F z = event.z(y, tan);
    F angle = compat::atan(tan);
    F cos = compat::cos(angle);

    F sec = 1 / cos;
    F sec_sq = sec * sec;

    F A = (((4 / (cos * cos)) * detector.inv_pow_sigma_dl) +
           (2 * tan * tan * detector.inv_pow_sigma_z));
    F B = -4 * tan * detector.inv_pow_sigma_z;
    F C = 2 * detector.inv_pow_sigma_z;
    F B_2 = (B / 2) * (B / 2);

    F bb_y = bby(A, C, B_2);
    F bb_z = bbz(A, C, B_2);

    Pixel<> center_pixel = detector.pixel_location(y, z);

    // bounding box limits for event
    Pixel<> ur(center_pixel.x - pixels_in_line(bb_y, detector.pixel_height),
               center_pixel.y + pixels_in_line(bb_z, detector.pixel_width));
    Pixel<> dl(center_pixel.x + pixels_in_line(bb_y, detector.pixel_height),
               center_pixel.y - pixels_in_line(bb_z, detector.pixel_width));
    Pixel<> ul(center_pixel.x - pixels_in_line(bb_y, detector.pixel_height),
               center_pixel.y - pixels_in_line(bb_z, detector.pixel_width));

    int bb_size = ceilf(((ur.y - ul.y) * (dl.x - ur.x)) / WARP_SIZE);

    int pixel_count = 0;

    Pixel<> first_pixel = ul;

    for (int k = 0; k < bb_size; ++k) {
      Pixel<> pixel = warp_space_pixel(offset, first_pixel, ul, ur, dl, tid);
      Point<F> point = detector.pixel_center(pixel);

      if (in_ellipse(A, B, C, y, z, point)) {
        point.x -= y;
        point.y -= z;

        F event_kernel = kernel(y,
                                tan,
                                sec,
                                sec_sq,
                                detector.radius,
                                point,
                                detector.inv_cor_mat_diag,
                                sqrt_det_cor_mat) /
                         tex2D<F>(sensitivity, pixel.x, pixel.y);

        acc += event_kernel * tex2D<F>(sensitivity, pixel.x, pixel.y) *
               rho[IMAGE_SPACE_LINEAR_INDEX(pixel.x, pixel.y)];

        sh_mem_pixel_buffer[SH_MEM_INDEX(threadIdx.x, pixel_count, 0)] =
            pixel.x;
        sh_mem_pixel_buffer[SH_MEM_INDEX(threadIdx.x, pixel_count, 1)] =
            pixel.y;
        pixel_count++;
      }
    }

    for (int xor_iter = 16; xor_iter >= 1; xor_iter /= 2) {
      acc += __shfl_xor(acc, xor_iter, WARP_SIZE);
    }

    F inv_acc = 1 / acc;

    for (int k = 0; k < pixel_count; ++k) {

      Pixel<> pixel(sh_mem_pixel_buffer[SH_MEM_INDEX(threadIdx.x, k, 0)],
                    sh_mem_pixel_buffer[SH_MEM_INDEX(threadIdx.x, k, 1)]);
      Point<F> point = detector.pixel_center(pixel);
      point.x -= y;
      point.y -= z;

      F event_kernel = kernel(y,
                              tan,
                              sec,
                              sec_sq,
                              detector.radius,
                              point,
                              detector.inv_cor_mat_diag,
                              sqrt_det_cor_mat) /
                       tex2D<F>(sensitivity, pixel.x, pixel.y);

      atomicAdd(
          &output_image[BUFFER_LINEAR_INDEX(pixel.x, pixel.y)],
          (event_kernel * rho[IMAGE_SPACE_LINEAR_INDEX(pixel.x, pixel.y)]) *
              inv_acc * sec_sq);
    }
  }
}
