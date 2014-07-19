#pragma once

#define WARP_SIZE 32

// reconstruction mode (comment out both for simple kernel)
#define THREAD_GRANULARITY 0  // single thread processes single event
#define WARP_GRANULARITY 1    // whole warp processes single event
#define USE_TEXTURE 1         // use regular CC 2.x compatible texture
#define USE_TEXTURE_OBJECT 0  // requires CC 3.0
#define SHARED_BUFFER 1       // shared memory pixel buffer in error ellipse
#define NORMAL_PHANTOM 0

#define PIXEL_INDEX(p) ((p.x * detector.n_z_pixels) + p.y)

#define COUNT_MAX_PIXELS_PER_THREAD 0  // count max pixels in thread
#define MAX_PIXELS_PER_THREAD 11       // this has been chosen based on above
#define MAX_THREADS_PER_BLOCK 512      // this has been chosen arbitrarily

#if USE_TEXTURE_OBJECT
#define TEX_ARG(v) cudaTextureObject_t v
#define TEX_VAL(v) v
#else
#define TEX_ARG(v) int v##_dummy
#define TEX_VAL(v) 0
#endif

#if USE_TEXTURE_OBJECT
#define TEX_2D(F, t, p) tex2D<F>(t, p.x, p.y)
#elif USE_TEXTURE
#define TEX_2D(F, t, p) tex2D(tex_##t, p.x, p.y)
#else
#define TEX_2D(F, t, p) 1
#endif
