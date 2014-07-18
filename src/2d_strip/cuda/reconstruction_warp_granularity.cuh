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

__device__ Pixel<> warp_space_pixel(int offset,
                                    Pixel<> ul,
                                    int width,
                                    int& index) __device__;
#ifdef SHARED_REGISTER
__shared__ float sqrt_det_cor_mat;
__shared__ int block_size;
__shared__ int number_of_blocks;
#endif

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

#ifdef SHARED_REGISTER
  if (blockIdx.x == 0) {

    sqrt_det_cor_mat = detector.sqrt_det_cor_mat();
    block_size = (n_blocks * (n_threads_per_block / WARP_SIZE));
    number_of_blocks = (n_events + block_size - 1) / block_size;
    __syncthreads();
  }
#else
  F sqrt_det_cor_mat = detector.sqrt_det_cor_mat();
  int block_size = (n_blocks * (n_threads_per_block / WARP_SIZE));
  int number_of_blocks = (n_events + block_size - 1) / block_size;

#endif

#ifdef SHARED_BUFFOR
  __shared__ short sh_mem_pixel_buffer[20 * 512];
#endif

  int thread_warp_index = threadIdx.x / WARP_SIZE;

  for (int i = 0; i < number_of_blocks; ++i) {

    // calculate offset in event memory location for warp
    int warp_id =
        (i * block_size + (blockIdx.x * (n_threads_per_block / WARP_SIZE)) +
         thread_warp_index);

    if (warp_id >= n_events)
      break;

    Event<F> event(
        events_z_u[warp_id], events_z_d[warp_id], events_dl[warp_id]);

    F acc = 0;

    F tan, y, z;
    event.transform(detector.radius, tan, y, z);
    F angle = std::atan(tan);

    Point<F> ellipse_center(y, z);

    F sec, sec_sq, A, B, C, bb_y, bb_z;
    detector.ellipse_bb(angle, tan, sec, sec_sq, A, B, C, bb_y, bb_z);

    Pixel<> center_pixel = detector.pixel_location(y, z);

    // bounding box limits for event
    Pixel<> tl(center_pixel.x - n_pixels_in_line(bb_y, detector.pixel_height),
               center_pixel.y - n_pixels_in_line(bb_z, detector.pixel_width));
    Pixel<> br(center_pixel.x + n_pixels_in_line(bb_y, detector.pixel_height),
               center_pixel.y + n_pixels_in_line(bb_z, detector.pixel_width));

    int bb_width = br.y - tl.y;
    int bb_height = br.x - tl.x;
    int bb_size = bb_width * bb_height;

    F inv_bb_width = F(1) / bb_width;

#ifdef SHARED_BUFFOR
    int pixel_count = 0;
#endif

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

        acc += event_kernel * TEX_2D(F, sensitivity, pixel) *
               rho[IMAGE_SPACE_LINEAR_INDEX(pixel)];
#ifdef SHARED_BUFFOR
        sh_mem_pixel_buffer[SH_MEM_INDEX(threadIdx.x, pixel_count, 0)] =
            pixel.x;
        sh_mem_pixel_buffer[SH_MEM_INDEX(threadIdx.x, pixel_count, 1)] =
            pixel.y;
        pixel_count++;
#endif
      }
    }

    for (int xor_iter = 16; xor_iter >= 1; xor_iter /= 2) {
      acc += __shfl_xor(acc, xor_iter, WARP_SIZE);
    }

    F inv_acc = 1 / acc;


#ifdef SHARED_BUFFOR
    for (int k = 0; k < pixel_count; ++k) {

      Pixel<> pixel(sh_mem_pixel_buffer[SH_MEM_INDEX(threadIdx.x, k, 0)],
                    sh_mem_pixel_buffer[SH_MEM_INDEX(threadIdx.x, k, 1)]);
      Point<F> point = detector.pixel_center(pixel);
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

      atomicAdd(&output[BUFFER_LINEAR_INDEX(pixel)],
                event_kernel * rho[IMAGE_SPACE_LINEAR_INDEX(pixel)] * inv_acc *
                    sec_sq);
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

        atomicAdd(
            &output[BUFFER_LINEAR_INDEX(pixel)],
            event_kernel * rho[IMAGE_SPACE_LINEAR_INDEX(pixel)] * inv_acc);
      }
    }
#endif
  }
}

__device__ Pixel<> warp_space_pixel(int offset,
                                    Pixel<> tl,
                                    int width,
                                    float inv_bb_width,
                                    int& index) {
  index = (threadIdx.x & (WARP_SIZE - 1)) + offset;
  Pixel<> pixel;
  pixel.x = index * inv_bb_width;
  pixel.y = index - width * pixel.x;
  pixel.x += tl.x;
  pixel.y += tl.y;
  return pixel;
}
