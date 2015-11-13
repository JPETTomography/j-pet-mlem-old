#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

#include "util/cuda/debug.h"  // catches all CUDA errors
#include "util/cuda/memory.h"
#include "../event.h"
#include "../kernel.h"
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
#define reconstruction reconstruction<Kernel><<<blocks, threads>>>
#else
  (void)n_blocks, n_threads_per_block;  // mark used
#define reconstruction reconstruction<Kernel>
#endif

  cudaSetDevice(device);
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device);
  device_name(prop.name);

  size_t image_size = scanner.total_n_pixels * sizeof(F);

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
  util::cuda::memory<F> output_rho(image_size);
  Output rho_output(width, height, output_rho.host_ptr);
  for (int i = 0; i < scanner.total_n_pixels; ++i) {
    output_rho[i] = 100;
  }
  output_rho.copy_to_device();

  // this class allocated CUDA pointers and deallocated them in destructor
  GPU::ResponsesSOA<F> responses_soa(responses, n_responses);

  for (int ib = 0; ib < n_iteration_blocks; ++ib) {
    for (int it = 0; it < n_iterations_in_block; ++it) {
      progress(ib * n_iterations_in_block + it, false);

      rho = output_rho;
      output_rho.zero_on_device();

      reconstruction(scanner,
                     responses_soa.z_u,
                     responses_soa.z_d,
                     responses_soa.dl,
                     n_responses,
                     output_rho.device_ptr);
      cudaThreadSynchronize();

      progress(ib * n_iterations_in_block + it, true);

      // always output first 5 iterations, and at 10, 15, 20, 30, 50, 100
      if (!ib && it < n_iterations_in_block - 1 &&
          (it < 5 || it == 9 || it == 14 || it == 19 || it == 29 || it == 49 ||
           it == 99)) {
        output_rho.copy_from_device();
        output(it + 1, rho_output);
      }
    }

    output_rho.copy_from_device();
    output((ib + 1) * n_iterations_in_block, rho_output);
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
