#pragma once

#define WARP_SIZE 32

// reconstruction mode (comment out both for simple kernel)
#define EVENT_GRANULARITY 1
//#define WARP_GRANULARITY 1

#define NORMAL_PHANTOM 0

#define IMAGE_SPACE_LINEAR_INDEX(Y, Z) (Y * detector.n_z_pixels) + Z
#define BUFFER_LINEAR_INDEX(Y, Z) \
  (blockIdx.x * detector.total_n_pixels) + (Y * detector.n_z_pixels) + Z
#define SH_MEM_INDEX(ID, N, I) (ID * 20 + (2 * N + I))
