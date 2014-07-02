#pragma once

#include <cuda_runtime.h>
#include <cmath>

#include "geometry/pixel.h"
#include "geometry/point.h"

#include "config.h"

#if OLD_WARP_SPACE_PIXEL
__device__ void warp_space_pixel(Pixel<>& pixel,
                                 int offset,
                                 int ul,
                                 int width,
                                 int height,
                                 int& index) {

  index = threadIdx.x & 31 + offset;
  pixel.y = index / width;
  pixel.x = index - width * pixel.y;
  pixel.y += ul.y;
  pixel.x += ul.x;
}
#else
__device__ Pixel<> warp_space_pixel(int& offset,
                                    Pixel<>& first_pixel,
                                    Pixel<> ul,
                                    Pixel<> ur,
                                    Pixel<> dl,
                                    int tid) {
  offset = ur.y - first_pixel.y;
  int tid_in_warp = threadIdx.x & 31;

  Pixel<> pixel;
  if (tid_in_warp < offset) {
    pixel = Pixel<>(first_pixel.x, first_pixel.y + tid_in_warp);
  } else {
    pixel = Pixel<>(
        first_pixel.x + (int)(__fdividef(tid - offset, ur.y - dl.y)) + 1,
        dl.y + ((tid - offset) % (ur.y - dl.y)));
  }

  first_pixel.y = ul.y + (32 - offset) % (ur.y - ul.y);
  first_pixel.x += int(__fdividef(32 - offset, ur.y - ul.y)) + 1;

  return pixel;
}
#endif

template <typename F> int n_pixels_in_line(F length, F pixel_size) $ {
  return int((length + 0.5f) / pixel_size);
}
