#pragma once

#include <cuda_runtime.h>

#include "2d/barrel/simple_geometry.h"
#include "common/types.h"

namespace Common {
namespace GPU {

using SimpleGeometry = PET2D::Barrel::SimpleGeometry<F, S, Hit>;
using Pixel = PET2D::Pixel<S>;

/// Calculates sensitivity out of given pixel_infos
__global__ void reduce_to_sensitivity(const Pixel* pixels,
                                      const F* pixel_weights,
                                      const size_t n_pixel_infos,
                                      float* output,
                                      const int width);

/// Invert given values i -> 1/i
__global__ void invert(float* input_output, const size_t size);

/// Reduce given register within the warp
template <typename Type> __device__ inline void reduce(Type& value) {
#if __CUDA_ARCH__ >= 300
  // reduce acc from all threads using __shfl_xor
  for (int i = 16; i >= 1; i /= 2) {
    value += __shfl_xor(value, i, WARP_SIZE);
  }
#else
  // fallback to older reduction algorithm
  __shared__ Type accumulator[MAX_THREADS_PER_BLOCK];
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
}  // Common
