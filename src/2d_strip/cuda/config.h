#pragma once

#define WARP_SIZE 32

// reconstruction mode (comment out both for simple kernel)
#define EVENT_GRANULARITY 1
//#define WARP_GRANULARITY 1
#define USE_TEXTURE 1
//#define USE_TEXTURE_OBJECT 1 // requires CC 3.0

#define NORMAL_PHANTOM 0

#define IMAGE_SPACE_LINEAR_INDEX(Y, Z) (Y * detector.n_z_pixels) + Z
#define BUFFER_LINEAR_INDEX(Y, Z) \
  (blockIdx.x * detector.total_n_pixels) + (Y * detector.n_z_pixels) + Z
#define SH_MEM_INDEX(ID, N, I) (ID * 20 + (2 * N + I))

#if USE_TEXTURE_OBJECT
#define TEX_ARG(v) cudaTextureObject_t v
#define TEX_VAL(v) v
#else
#define TEX_ARG(v) int v##_dummy
#define TEX_VAL(v) 0
#endif

#if USE_TEXTURE_OBJECT
#define TEX_2D(F, t, y, z) tex2D<F>(t, y, z)
#elif USE_TEXTURE
#define TEX_2D(F, t, y, z) tex2D(tex_##t, y, z)
#else
#define TEX_2D(F, t, y, z) 1
#endif
