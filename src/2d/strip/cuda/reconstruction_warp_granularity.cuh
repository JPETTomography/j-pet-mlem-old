#pragma once

#include <cuda_runtime.h>

#include "2d/geometry/point.h"
#include "../event.h"
#include "../detector.h"

#define PIXEL_INDEX(p) (((p).y * detector.n_z_pixels) + (p).x)

namespace PET2D {
namespace Strip {
namespace GPU {

template <typename F> __device__ void reduce(F& value);

template <template <typename Float> class Kernel, typename F>
__global__ void reconstruction(Detector<F> detector,
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

  const int n_warps_per_block = n_threads_per_block / WARP_SIZE;
  const int n_warps = n_blocks * n_warps_per_block;
  const int max_events_per_warp = (n_events + n_warps - 1) / n_warps;
  const int warp_index = threadIdx.x / WARP_SIZE;

#if CACHE_ELLIPSE_PIXELS
  // gathers all pixel coordinates inside 3 sigma ellipse
  __shared__ short2
      ellipse_pixels[MAX_PIXELS_PER_THREAD][MAX_THREADS_PER_BLOCK];
#endif

  Kernel<F> kernel(detector.sigma_z, detector.sigma_dl);

  for (int i = 0; i < max_events_per_warp; ++i) {

    int event_index = i * n_warps + blockIdx.x * n_warps_per_block + warp_index;

    if (event_index >= n_events)
      break;

    Event event(events_z_u[event_index],
                events_z_d[event_index],
                events_dl[event_index]);

    F denominator = 0;

    F tan, y, z;
    event.transform(detector.radius, tan, y, z);

    F sec, A, B, C, bb_y, bb_z;
    kernel.ellipse_bb(tan, sec, A, B, C, bb_y, bb_z);

    Point ellipse_center(z, y);
    Pixel center_pixel = detector.pixel_at(ellipse_center);

    // bounding box limits for event
    const int bb_half_width = n_pixels_in_line(bb_z, detector.pixel_width);
    const int bb_half_height = n_pixels_in_line(bb_y, detector.pixel_height);
    Pixel bb_tl(center_pixel.x - bb_half_width,
                center_pixel.y - bb_half_height);
    Pixel bb_br(center_pixel.x + bb_half_width,
                center_pixel.y + bb_half_height);
    Pixel detector_tl(0, 0);
    Pixel detector_br(detector.n_z_pixels - 1, detector.n_y_pixels - 1);

    // check boundary conditions
    bb_tl.clamp(detector_tl, detector_br);
    bb_br.clamp(detector_tl, detector_br);

    const int bb_width = bb_br.x - bb_tl.x;
    const int bb_height = bb_br.y - bb_tl.y;
    const int bb_size = bb_width * bb_height;
    F inv_bb_width = F(1) / bb_width;

#if CACHE_ELLIPSE_PIXELS
    int n_ellipse_pixels = 0;
    F ellipse_kernel_mul_rho[MAX_PIXELS_PER_THREAD];
#endif

    for (int offset = 0; offset < bb_size; offset += WARP_SIZE) {
      int index;
      Pixel pixel =
          warp_space_pixel(offset, bb_tl, bb_width, inv_bb_width, index);

      if (index >= bb_size)
        break;

      Point point = detector.pixel_center(pixel);

      if (kernel.in_ellipse(A, B, C, ellipse_center, point)) {
        point -= ellipse_center;

        F pixel_sensitivity =
            USE_SENSITIVITY ? tex2D(tex_sensitivity, pixel.x, pixel.y) : 1;

        F event_kernel =
            USE_KERNEL ? kernel(y, tan, sec, detector.radius, point) : 1;

        F event_kernel_mul_rho =
            event_kernel * tex2D(tex_rho, pixel.x, pixel.y);
        denominator += event_kernel_mul_rho * pixel_sensitivity;

#if CACHE_ELLIPSE_PIXELS
        ellipse_pixels[n_ellipse_pixels][threadIdx.x] =
            make_short2(pixel.x, pixel.y);
        ellipse_kernel_mul_rho[n_ellipse_pixels] = event_kernel_mul_rho;
        ++n_ellipse_pixels;
#endif
      }
    }

    // reduce denominator so all threads now share same value
    reduce(denominator);

    F inv_acc = 1 / denominator;

#if CACHE_ELLIPSE_PIXELS
    for (int p = 0; p < n_ellipse_pixels; ++p) {
      short2 pixel = ellipse_pixels[p][threadIdx.x];

      atomicAdd(&output_rho[PIXEL_INDEX(pixel)],
                ellipse_kernel_mul_rho[p] * inv_acc);
    }
#else
    for (int offset = 0; offset < bb_size; offset += WARP_SIZE) {
      int index;
      Pixel pixel =
          warp_space_pixel(offset, top_left, bb_width, inv_bb_width, index);

      if (index >= bb_size)
        break;

      Point point = detector.pixel_center(pixel);

      if (detector.in_ellipse(A, B, C, ellipse_center, point)) {
        point -= ellipse_center;

        F event_kernel =
            USE_KERNEL ? kernel(y, tan, sec, detector.radius, point) : 1;

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

template <typename F> __device__ void reduce(F& value) {
#if __CUDA_ARCH__ >= 300
  // reduce acc from all threads using __shfl_xor
  for (int i = 16; i >= 1; i /= 2) {
    value += __shfl_xor(value, i, WARP_SIZE);
  }
#else
  // fallback to older reduction algorithm
  __shared__ F accumulator[MAX_THREADS_PER_BLOCK];
  int tid = threadIdx.x;
  int index = (tid & (WARP_SIZE - 1));
  accumulator[tid] = value;
  for (int i = 16; i >= 1; i /= 2) {
    if (index < i)
      accumulator[tid] += accumulator[tid + i];
  }
  value = accumulator[tid & ~(WARP_SIZE - 1)];
#endif
}
}  // GPU
}  // Strip
}  // PET2D
