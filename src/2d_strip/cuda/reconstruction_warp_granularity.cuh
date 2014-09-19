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

#define ALWAYS_TRUE /***/ (n_blocks > 0)
#define ALWAYS_FALSE /**/ (n_blocks == 0)

  // set conditions that always evaluate to same value, if we use kernel
  // parameter that always evaluate to same value, compiler cannot optimize it
  // out and must generate both branches
  bool use_kernel = true;       // use ALWAYS_XXX to measure performance
  bool use_sensitivity = true;  // use ALWAYS_XXX to measure performance

  float full_acc = 0;

  const F sqrt_det_cor_mat = detector.sqrt_det_cor_mat();
  const int n_warps_per_block = n_threads_per_block / WARP_SIZE;
  const int n_warps = n_blocks * n_warps_per_block;
  const int max_events_per_warp = (n_events + n_warps - 1) / n_warps;
  const int warp_index = threadIdx.x / WARP_SIZE;

#if CACHE_ELLIPSE_PIXELS
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

    Point<F> ellipse_center(z, y);

    F sec, sec_sq, A, B, C, bb_y, bb_z;
    detector.ellipse_bb(angle, tan, sec, sec_sq, A, B, C, bb_y, bb_z);

    Pixel<> center_pixel = detector.pixel_location(ellipse_center);

    // bounding box limits for event
    const int bb_half_width = n_pixels_in_line(bb_z, detector.pixel_width);
    const int bb_half_height = n_pixels_in_line(bb_y, detector.pixel_height);
    const Pixel<> top_left(center_pixel.x - bb_half_width,
                           center_pixel.y - bb_half_height);

    const int bb_width = 2 * bb_half_width;
    const int bb_height = 2 * bb_half_height;
    const int bb_size = bb_width * bb_height;
    F inv_bb_width = F(1) / bb_width;

#if CACHE_ELLIPSE_PIXELS
    int n_ellipse_pixels = 0;
    F ellipse_kernel_mul_rho[MAX_PIXELS_PER_THREAD];
#endif

    for (int offset = 0; offset < bb_size; offset += WARP_SIZE) {
      int index;
      Pixel<> pixel =
          warp_space_pixel(offset, top_left, bb_width, inv_bb_width, index);

      if (index >= bb_size)
        break;

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

        F event_kernel_mul_rho =
            event_kernel * tex2D(tex_rho, pixel.x, pixel.y);
        acc += event_kernel_mul_rho * pixel_sensitivity;
#if CACHE_ELLIPSE_PIXELS
        ellipse_pixels[n_ellipse_pixels][threadIdx.x] =
            make_short2(pixel.x, pixel.y);
        ellipse_kernel_mul_rho[n_ellipse_pixels] = event_kernel_mul_rho;
        ++n_ellipse_pixels;
#endif
      }
    }

#if __CUDA_ARCH__ >= 300
    // reduce acc from all threads using __shfl_xor
    for (int xor_iter = 16; xor_iter >= 1; xor_iter /= 2) {
      acc += __shfl_xor(acc, xor_iter, WARP_SIZE);
    }
#else
    // fallback to older reduction algorithm
    __shared__ F accumulator[MAX_THREADS_PER_BLOCK];
    int tid = threadIdx.x;
    int index = (tid & (WARP_SIZE - 1));
    accumulator[tid] = acc;
    if (index < 16)
      accumulator[tid] += accumulator[tid + 16];
    if (index < 8)
      accumulator[tid] += accumulator[tid + 8];
    if (index < 4)
      accumulator[tid] += accumulator[tid + 4];
    if (index < 2)
      accumulator[tid] += accumulator[tid + 2];
    if (index < 1)
      accumulator[tid] += accumulator[tid + 1];
    acc = accumulator[tid & ~(WARP_SIZE - 1)];
#endif

    F inv_acc = 1 / acc;
    full_acc += acc;

#if CACHE_ELLIPSE_PIXELS
    for (int p = 0; p < n_ellipse_pixels; ++p) {
      short2 pixel = ellipse_pixels[p][threadIdx.x];

#if !NO_ATOMIC
        atomicAdd(&output_rho[PIXEL_INDEX(pixel)],
                  ellipse_kernel_mul_rho[p] * inv_acc);
#else
        output_rho[PIXEL_INDEX(pixel)] += ellipse_kernel_mul_rho[p] * inv_acc;
#endif
    }
#else
    for (int offset = 0; offset < bb_size; offset += WARP_SIZE) {
      int index;
      Pixel<> pixel =
          warp_space_pixel(offset, top_left, bb_width, inv_bb_width, index);

      if (index >= bb_size)
        break;

      Point<F> point = detector.pixel_center(pixel);

      if (detector.in_ellipse(A, B, C, ellipse_center, point)) {
        point -= ellipse_center;

#if USE_SENSITIVITY
        F inv_pixel_sensitivity = tex2D(tex_inv_sensitivity, pixel.x, pixel.y);
#else
        F inv_pixel_sensitivity = 1;
#endif

        F event_kernel = kernel(y,
                                tan,
                                sec,
                                sec_sq,
                                detector.radius,
                                point,
                                detector.inv_cor_mat_diag,
                                sqrt_det_cor_mat) *
                         inv_pixel_sensitivity;

        atomicAdd(&output_rho[PIXEL_INDEX(pixel)],
                  event_kernel * tex2D(tex_rho, pixel.x, pixel.y) * inv_acc);
      }
    }
#endif
  }
}

template <typename F> _ int n_pixels_in_line(F length, F pixel_size) {
  return (length + F(0.5)) / pixel_size;
}

template <typename F>
__device__ Pixel<> warp_space_pixel(int offset,
                                    Pixel<> tl,
                                    int width,
                                    F inv_width,
                                    int& index) {
  // threadIdx.x % WARP_SIZE + offset : works for WARP_SIZE = 2^n
  index = (threadIdx.x & (WARP_SIZE - 1)) + offset;
  Pixel<> pixel;
  pixel.y = index * inv_width;  // index/width but faster
  pixel.x = index - width * pixel.y;
  pixel.x += tl.x;
  pixel.y += tl.y;
  return pixel;
}
