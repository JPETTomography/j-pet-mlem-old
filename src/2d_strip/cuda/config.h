#pragma once

#define WARP_SIZE 32

// reconstruction mode (comment out both for simple kernel)
#define EVENT_GRANULARITY 1
//#define WARP_GRANULARITY 1
#define USE_TEXTURE 1
//#define USE_TEXTURE_OBJECT 1 // requires CC 3.0
#define SHARED_REGISTER 1  // nearly 4x speedup
//#define SHARED_BUFFER 1 // shared memory pixel buffer in error ellipse

#define NORMAL_PHANTOM 0

#define IMAGE_SPACE_LINEAR_INDEX(p) (p.x * detector.n_z_pixels) + p.y
#define BUFFER_LINEAR_INDEX(p) \
  (blockIdx.x * detector.total_n_pixels) + IMAGE_SPACE_LINEAR_INDEX(p)
#define SH_MEM_INDEX(ID, N, I) (ID * 20 + (2 * N + I))

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
