#include "kernels.h"

namespace Common {
namespace GPU {

__global__ void reduce_to_sensitivity(const Pixel* pixels,
                                      const F* pixel_weights,
                                      const size_t n_pixel_infos,
                                      float* output,
                                      const int width) {

  const auto tid = (blockIdx.x * blockDim.x) + threadIdx.x;
  const auto n_threads = gridDim.x * blockDim.x;
  const auto n_chunks = (n_pixel_infos + n_threads - 1) / n_threads;

  for (int chunk = 0; chunk < n_chunks; ++chunk) {
    int index = chunk * n_threads + tid;
    // check if we are still on the list
    if (index >= n_pixel_infos) {
      break;
    }
    const auto pixel = pixels[index];
    atomicAdd(&output[pixel.y * width + pixel.x], pixel_weights[index]);
  }
}

__global__ void invert(float* input_output, const size_t size) {

  const auto tid = (blockIdx.x * blockDim.x) + threadIdx.x;
  const auto n_threads = gridDim.x * blockDim.x;
  const auto n_chunks = (size + n_threads - 1) / n_threads;

  for (int chunk = 0; chunk < n_chunks; ++chunk) {
    int index = chunk * n_threads + tid;
    // check if we are still on the list
    if (index >= size) {
      break;
    }
    auto input = input_output[index];
    if (input > 0) {
      input_output[index] = 1 / input;
    }
  }
}

}  // GPU
}  // Common
