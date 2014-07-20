#pragma once

#include <cuda_runtime.h>

#include "geometry/point.h"
#include "../event.h"
#include "../kernel.h"
#include "../strip_detector.h"

#include "config.h"

template <typename F>
__global__ void reconstruction(StripDetector<F> detector,
                               F* events_z_u,
                               F* events_z_d,
                               F* events_dl,
                               const int n_events,
                               F* output_rho,
                               const int n_blocks,
                               const int n_threads_per_block) {
  Kernel<F> kernel;

  const F sqrt_det_cor_mat = detector.sqrt_det_cor_mat();
  const int n_warps_per_block = n_threads_per_block / WARP_SIZE;
  const int n_warps = n_blocks * n_warps_per_block;
  const int max_events_per_warp = (n_events + n_warps - 1) / n_warps;
  const int warp_index = threadIdx.x / WARP_SIZE;

#if SHARED_BUFFER
  // gathers all pixel coordinates inside 3 sigma ellipse
  __shared__ short2
      ellipse_pixels[MAX_PIXELS_PER_THREAD][MAX_THREADS_PER_BLOCK];
#endif

  for (int i = 0; i < max_events_per_warp; ++i) {

    int event_index = i * n_warps + blockIdx.x * n_warps_per_block + warp_index;

    if (event_index >= n_events)
      break;

    Event<F> event(events_z_u[event_index],
                   events_z_d[event_index],
                   events_dl[event_index]);

    F acc = 0;

    F tan, y, z;
    event.transform(detector.radius, tan, y, z);
    F angle = std::atan(tan);

    Point<F> ellipse_center(y, z);

    F sec, sec_sq, A, B, C, bb_y, bb_z;
    detector.ellipse_bb(angle, tan, sec, sec_sq, A, B, C, bb_y, bb_z);

    Pixel<> center_pixel = detector.pixel_location(y, z);

    // bounding box limits for event
    const int bb_half_width = n_pixels_in_line(bb_z, detector.pixel_width);
    const int bb_half_height = n_pixels_in_line(bb_y, detector.pixel_height);
    const Pixel<> tl(center_pixel.x - bb_half_height,
                     center_pixel.y - bb_half_width);

    const int bb_width = 2 * bb_half_width;
    const int bb_height = 2 * bb_half_height;
    const int bb_size = bb_width * bb_height;
    F inv_bb_width = F(1) / bb_width;

#if SHARED_BUFFER
    int n_ellipse_pixels = 0;
#endif

    F ellipse_kernel_mul_rho[MAX_PIXELS_PER_THREAD];

    for (int offset = 0; offset < bb_size; offset += WARP_SIZE) {
      int index;
      Pixel<> pixel =
          warp_space_pixel(offset, tl, bb_width, inv_bb_width, index);

      if (index >= bb_size)
        break;

      Point<F> point = detector.pixel_center(pixel);

      if (detector.in_ellipse(A, B, C, ellipse_center, point)) {
        point -= ellipse_center;

#if SENSITIVITY_TEXTURE
        F pixel_sensitivity = tex2D(tex_sensitivity, pixel.y, pixel.x);
#else
        F pixel_sensitivity = 1;
#endif
        F event_kernel = kernel(y,
                                tan,
                                sec,
                                sec_sq,
                                detector.radius,
                                point,
                                detector.inv_cor_mat_diag,
                                sqrt_det_cor_mat) /
                         pixel_sensitivity;

        F event_kernel_mul_rho =
            event_kernel * tex2D(tex_rho, pixel.y, pixel.x);
        acc += event_kernel_mul_rho * pixel_sensitivity;
#if SHARED_BUFFER
        ellipse_pixels[n_ellipse_pixels][threadIdx.x] =
            make_short2(pixel.x, pixel.y);
        ellipse_kernel_mul_rho[n_ellipse_pixels] = event_kernel_mul_rho;
        ++n_ellipse_pixels;
#endif
      }
    }

    for (int xor_iter = 16; xor_iter >= 1; xor_iter /= 2) {
      acc += __shfl_xor(acc, xor_iter, WARP_SIZE);
    }

    F inv_acc = 1 / acc;

#if SHARED_BUFFER
    for (int p = 0; p < n_ellipse_pixels; ++p) {
      short2 pixel = ellipse_pixels[p][threadIdx.x];
      F event_kernel_mul_rho = ellipse_kernel_mul_rho[p];
      atomicAdd(&output_rho[PIXEL_INDEX(pixel)],
                event_kernel_mul_rho * inv_acc * sec_sq);
    }
#else
    for (int offset = 0; offset < bb_size; offset += WARP_SIZE) {
      int index;
      Pixel<> pixel =
          warp_space_pixel(offset, tl, bb_width, inv_bb_width, index);

      if (index >= bb_size)
        break;

      Point<F> point = detector.pixel_center(pixel);

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
                         TEX_2D(F, sensitivity, pixel);

        atomicAdd(&output_rho[PIXEL_INDEX(pixel)],
                  event_kernel * rho[PIXEL_INDEX(pixel)] * inv_acc);
      }
    }
#endif
  }
}

template <typename F>
__device__ Pixel<> warp_space_pixel(int offset,
                                    Pixel<> tl,
                                    int width,
                                    F inv_width,
                                    int& index) {
  index = (threadIdx.x & (WARP_SIZE - 1)) + offset;
  Pixel<> pixel;
  pixel.x = index * inv_width;
  pixel.y = index - width * pixel.x;
  pixel.x += tl.x;
  pixel.y += tl.y;
  return pixel;
}

template <typename F> _ int n_pixels_in_line(F length, F pixel_size) {
  return (length + F(0.5)) / pixel_size;
}
