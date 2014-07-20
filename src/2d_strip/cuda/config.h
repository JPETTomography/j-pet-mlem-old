#pragma once

#define WARP_SIZE 32

// reconstruction mode (comment out both for simple kernel)
#define THREAD_GRANULARITY 0    // single thread processes single event
#define WARP_GRANULARITY 1      // whole warp processes single event
#define USE_SENSITIVITY 1       // use sensitivity texture
#define CACHE_ELLIPSE_PIXELS 1  // shared memory pixel buffer in error ellipse
#define NORMAL_PHANTOM 0

// warp granulaty specific
#define MAX_PIXELS_PER_THREAD 12   // this has been chosen arbitrarily,
                                   // however in current implementation we have
                                   // no more than 11 pixels per thread.
#define MAX_THREADS_PER_BLOCK 512  // more does not make sense anyway

// generic macro returning pixel location in linear memory
#define PIXEL_INDEX(p) (((p).y * detector.n_z_pixels) + (p).x)
