#pragma once

#define WARP_SIZE 32

// warp granulaty specific
#define MAX_PIXELS_PER_THREAD 12   // this has been chosen arbitrarily,
                                   // however in current implementation we have
                                   // no more than 11 pixels per thread.
#define MAX_THREADS_PER_BLOCK 512  // more does not make sense anyway

// generic macro returning pixel location in linear memory

#define PIXEL_INDEX(p) (((p).y * detector.n_z_pixels) + (p).x)
#define WARP_BUFFER_PIXEL_INDEX(p) \
  (warp_stride + ((p).y * detector.n_z_pixels) + (p).x)
