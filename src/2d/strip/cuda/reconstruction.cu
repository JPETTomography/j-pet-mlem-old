#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

#include "util/cuda/debug.h"  // catches all CUDA errors
#include "util/cuda/memory.h"
#include "../event.h"
#include "../analytic_kernel.h"
#include "gpu_events_soa.h"
#include "reconstruction.h"

texture<float, 2, cudaReadModeElementType> tex_sensitivity;
texture<float, 2, cudaReadModeElementType> tex_rho;

#if USE_WARP_GRANULARITY
#include "reconstruction/warp_granularity.cuh"
#elif USE_THREAD_GRANULARITY
#include "reconstruction/thread_granularity.cuh"
#else
#include "reconstruction/simple.cuh"
#endif

namespace PET2D {
namespace Strip {
namespace GPU {
namespace Reconstruction {

template <typename F>
void fill_with_sensitivity(F* sensitivity, Scanner<F, S>& scanner);

template <typename F>
void run(Scanner<F, S>& scanner,
         Response<F>* responses,
         int n_responses,
         Output& rho_output,
         int n_iteration_blocks,
         int n_iterations_in_block,
         util::delegate<void(int iteration, const Output& output)> output,
         util::delegate<void(int completed, bool finished)> progress,
         int device,
         int n_blocks,
         int n_threads_per_block,
         util::delegate<void(const char* device_name)> device_name) {

#if __CUDACC__
  dim3 blocks(n_blocks);
  dim3 threads(n_threads_per_block);
#define reconstruction reconstruction<AnalyticKernel><<<blocks, threads>>>
#else
  (void)n_blocks, n_threads_per_block;  // mark used
#define reconstruction reconstruction<AnalyticKernel>
#endif

  cudaSetDevice(device);
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device);
  device_name(prop.name);

  const size_t image_size = scanner.total_n_pixels;
  const int width = scanner.n_z_pixels;
  const int height = scanner.n_y_pixels;

  util::cuda::texture2D<F> sensitivity(tex_sensitivity, width, height);
  {
    F output_sensitivity[width * height];
    fill_with_sensitivity(output_sensitivity, scanner);
    sensitivity = output_sensitivity;
    Output sensitivity_output(width, height, output_sensitivity);
    output(-1, sensitivity_output);
  }

  util::cuda::texture2D<F> rho(tex_rho, width, height);
  util::cuda::memory<F> output_rho(image_size, rho_output.data);
  output_rho.copy_to_device();

  // this class allocated CUDA pointers and deallocated them in destructor
  GPU::ResponsesSOA<F> responses_soa(responses, n_responses);

  for (int block = 0; block < n_iteration_blocks; ++block) {
    for (int i = 0; i < n_iterations_in_block; ++i) {
      int iteration = block * n_iterations_in_block + i;
      progress(iteration, false);

      rho = output_rho;
      output_rho.zero_on_device();

      reconstruction(scanner,
                     responses_soa.z_u,
                     responses_soa.z_d,
                     responses_soa.dl,
                     n_responses,
                     output_rho.device_ptr);
      cudaThreadSynchronize();

      progress(iteration, true);

      // always output first 5 iterations, and at 10, 15, 20, 30, 50, 100
      if (!block && i < n_iterations_in_block - 1 &&
          (i < 5 || i == 9 || i == 14 || i == 19 || i == 29 || i == 49 ||
           i == 99)) {
        output_rho.copy_from_device();
        output(i + 1, rho_output);
      }
    }

    output_rho.copy_from_device();
    int iteration = (block + 1) * n_iterations_in_block;
    output(iteration, rho_output);
  }

  progress(n_iteration_blocks * n_iterations_in_block, false);
}

template <typename F>
void fill_with_sensitivity(F* sensitivity, Scanner<F, S>& scanner) {

  size_t width = scanner.n_z_pixels;
  size_t height = scanner.n_y_pixels;

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      sensitivity[y * width + x] = scanner.pixel_sensitivity(Pixel(x, y));
    }
  }
}

// Explicit template instantiation
template void run(
    Scanner<float, S>& scanner,
    Response<float>* responses,
    int n_responses,
    Output& rho_output,
    int n_iteration_blocks,
    int n_iterations_in_block,
    util::delegate<void(int iteration, const Output& output)> output,
    util::delegate<void(int completed, bool finished)> progress,
    int device,
    int n_blocks,
    int n_threads_per_block,
    util::delegate<void(const char* device_name)> device_name);

}  // Reconstruction
}  // GPU
}  // Strip
}  // PET2D
