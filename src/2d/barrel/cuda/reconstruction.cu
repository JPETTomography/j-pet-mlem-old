#include <cuda_runtime.h>

#include "util/cuda/debug.h"  // catches all CUDA errors
#include "util/cuda/memory.h"
#include "util/delegate.h"

#include "reconstruction.h"

namespace PET2D {
namespace Barrel {
namespace GPU {
namespace Reconstruction {

texture<float, 2, cudaReadModeElementType> tex_rho;

// foreach p: count y[p] and store it in output_rho[p]
__global__ static void kernel_1(const PixelInfo* pixel_infos,
                                const size_t* lor_pixel_info_start,
                                const size_t* lor_pixel_info_end,
                                const Mean* means,
                                const int n_means,
                                float* output_rho,
                                const int width,
                                const int n_blocks,
                                const int n_threads_per_block) {

  const auto tid = (blockIdx.x * blockDim.x) + threadIdx.x;

  const auto n_threads = n_blocks * n_threads_per_block;
  const auto n_chunks = (n_means + n_threads - 1) / n_threads;

  for (int chunk = 0; chunk < n_chunks; ++chunk) {
    int mean_index = chunk * n_threads + tid;
    // check if we are still on the list
    if (mean_index >= n_means) {
      break;
    }

    auto mean = means[mean_index];
    auto lor_index = mean.lor.index();
    auto pixel_info_start = lor_pixel_info_start[lor_index];
    auto pixel_info_end = lor_pixel_info_end[lor_index];

    // count u for current lor
    F u = 0;
    for (auto i = pixel_info_start; i < pixel_info_end; ++i) {
      auto pixel_info = pixel_infos[i];
      auto pixel = pixel_info.pixel;
      u += tex2D(tex_rho, pixel.x, pixel.y) * pixel_info.weight;
    }
    F phi = mean.mean / u;
    for (auto i = pixel_info_start; i < pixel_info_end; ++i) {
      auto pixel_info = pixel_infos[i];
      auto pixel = pixel_info.pixel;
      atomicAdd(&output_rho[pixel.y * width + pixel.x],
                phi * pixel_info.weight);
    }
  }
}

// foreach p: count output_rho[p] *= rho[p]
__global__ static void kernel_2(float* output_rho,
                                const float* scale,
                                const int width,
                                const int height,
                                const int n_blocks,
                                const int n_threads_per_block) {

  const auto tid = (blockIdx.x * blockDim.x) + threadIdx.x;

  const auto n_threads = n_blocks * n_threads_per_block;
  const auto n_pixels = width * height;
  const auto n_pixel_chunks = (n_pixels + n_threads - 1) / n_threads;

  for (int chunk = 0; chunk < n_pixel_chunks; ++chunk) {
    int pixel_index = chunk * n_threads + tid;
    // check if we are still on the list
    if (pixel_index >= n_pixels) {
      break;
    }
    Pixel pixel(pixel_index % width, pixel_index / width);
    // there is no collision there, so we don't need atomics
    output_rho[pixel_index] *=
        tex2D(tex_rho, pixel.x, pixel.y) * scale[pixel_index];
  }
}

// calculates sensitivity out of given pixel_infos
__global__ static void sensitivity_kernel(const PixelInfo* pixel_infos,
                                          const size_t n_pixel_infos,
                                          float* output,
                                          const int width,
                                          const int n_blocks,
                                          const int n_threads_per_block) {

  const auto tid = (blockIdx.x * blockDim.x) + threadIdx.x;

  const auto n_threads = n_blocks * n_threads_per_block;
  const auto n_chunks = (n_pixel_infos + n_threads - 1) / n_threads;

  for (int chunk = 0; chunk < n_chunks; ++chunk) {
    int index = chunk * n_threads + tid;
    // check if we are still on the list
    if (index >= n_pixel_infos) {
      break;
    }
    PixelInfo pixel_info = pixel_infos[index];
    auto pixel = pixel_info.pixel;
    atomicAdd(&output[pixel.y * width + pixel.x], pixel_info.weight);
  }
}

// inverts values in input_output
__global__ static void invert_kernel(float* input_output,
                                     const size_t size,
                                     const int n_blocks,
                                     const int n_threads_per_block) {

  const auto tid = (blockIdx.x * blockDim.x) + threadIdx.x;

  const auto n_threads = n_blocks * n_threads_per_block;
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

void run(const SimpleGeometry& geometry,
         const Mean* means,
         int n_means,
         int width,
         int height,
         int n_iteration_blocks,
         int n_iterations_in_block,
         util::delegate<void(int iteration, float* rho)> output,
         util::delegate<void(int completed, bool finished)> progress,
         int device,
         int n_blocks,
         int n_threads_per_block,
         util::delegate<void(const char* device_name)> info) {

#if __CUDACC__
  dim3 blocks(n_blocks);
  dim3 threads(n_threads_per_block);
#endif

  cudaSetDevice(device);
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device);
  info(prop.name);

  util::cuda::on_device<PixelInfo> device_pixel_infos(geometry.pixel_infos,
                                                      geometry.n_pixel_infos);
  util::cuda::on_device<size_t> device_lor_pixel_info_start(
      geometry.lor_pixel_info_start, geometry.n_lors);
  util::cuda::on_device<size_t> device_lor_pixel_info_end(
      geometry.lor_pixel_info_end, geometry.n_lors);
  util::cuda::on_device<Mean> device_means(means, n_means);

  util::cuda::memory2D<F> rho(tex_rho, width, height);
  util::cuda::on_device<F> output_rho((size_t)width * height);

  for (auto& v : rho) {
    v = 1;
  }
  rho.copy_to_device();

  util::cuda::on_device<F> scale((size_t)width * height);
#if __CUDACC__
#define sensitivity_kernel sensitivity_kernel<<<blocks, threads>>>
#define invert_kernel invert_kernel<<<blocks, threads>>>
#endif
  scale.zero_on_device();
  sensitivity_kernel(device_pixel_infos,
                     geometry.n_pixel_infos,
                     scale,
                     width,
                     n_blocks,
                     n_threads_per_block);
  cudaThreadSynchronize();
  invert_kernel(scale, width * height, n_blocks, n_threads_per_block);
  cudaThreadSynchronize();

  for (int ib = 0; ib < n_iteration_blocks; ++ib) {
    for (int it = 0; it < n_iterations_in_block; ++it) {
      progress(ib * n_iterations_in_block + it, false);

      output_rho.zero_on_device();

#if __CUDACC__
#define kernel_1 kernel_1<<<blocks, threads>>>
#define kernel_2 kernel_2<<<blocks, threads>>>
#endif

      kernel_1(device_pixel_infos,
               device_lor_pixel_info_start,
               device_lor_pixel_info_end,
               device_means,
               n_means,
               output_rho,
               width,
               n_blocks,
               n_threads_per_block);

      cudaThreadSynchronize();

      kernel_2(output_rho, scale, width, height, n_blocks, n_threads_per_block);

      cudaThreadSynchronize();

      rho = output_rho;
      progress(ib * n_iterations_in_block + it, true);
    }

    rho.copy_from_device();
    output((ib + 1) * n_iterations_in_block, rho.host_ptr);
  }

  progress(n_iteration_blocks * n_iterations_in_block, false);
}

}  // Reconstruction
}  // GPU
}  // Barrel
}  // PET2D
