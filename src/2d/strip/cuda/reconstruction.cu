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
#include "reconstruction_warp_granularity.cuh"
#elif USE_THREAD_GRANULARITY
#include "reconstruction_thread_granularity.cuh"
#else
#include "reconstruction_simple.cuh"
#endif

namespace PET2D {
namespace Strip {
namespace GPU {
namespace Reconstruction {

template <typename F>
void fill_with_sensitivity(F* sensitivity, Scanner<F, short>& scanner);

template <typename F>
void run(Scanner<F, short>& scanner,
         Response<F>* responses,
         int n_responses,
         int n_iteration_blocks,
         int n_iterations_in_block,
         util::delegate<void(int, F*)> output,
         util::delegate<void(int, bool)> progress,
         int device,
         int n_blocks,
         int n_threads_per_block,
         util::delegate<void(const char*)> device_name) {

  cudaSetDevice(device);
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device);
  device_name(prop.name);

#if __CUDACC__
  dim3 blocks(n_blocks);
  dim3 threads(n_threads_per_block);
#endif

  size_t image_size = scanner.total_n_pixels * sizeof(F);

  const int width = scanner.n_z_pixels;
  const int height = scanner.n_y_pixels;

  util::cuda::memory2D<F> sensitivity(tex_sensitivity, width, height);
  fill_with_sensitivity(sensitivity.host_ptr, scanner);
  output(-1, sensitivity.host_ptr);

  util::cuda::memory2D<F> rho(tex_rho, width, height);
  for (int i = 0; i < scanner.total_n_pixels; ++i) {
    rho[i] = 100;
  }

  // this class allocated CUDA pointers and deallocated them in destructor
  GPU::ResponsesSOA<F> responses_soa(responses, n_responses);

#if USE_RHO_PER_WARP
  util::cuda::memory<F> output_rho(n_blocks, image_size);
#else
  util::cuda::on_device<F> output_rho(image_size);
#endif

  for (int ib = 0; ib < n_iteration_blocks; ++ib) {
    for (int it = 0; it < n_iterations_in_block; ++it) {
      progress(ib * n_iterations_in_block + it, false);

      output_rho.zero_on_device();
      rho.copy_to_device();

#if __CUDACC__
#define reconstruction reconstruction<Kernel><<<blocks, threads>>>
#endif
      reconstruction(scanner,
                     responses_soa.z_u,
                     responses_soa.z_d,
                     responses_soa.dl,
                     n_responses,
                     output_rho.device_ptr,
                     n_blocks,
                     n_threads_per_block);

      cudaThreadSynchronize();

#if USE_RHO_PER_WARP
      output_rho.copy_from_device();

      for (int i = 0; i < scanner.n_y_pixels; ++i) {
        for (int j = 0; j < scanner.n_z_pixels; ++j) {
          int pixel_adr = i * scanner.n_y_pixels + j;
          rho[pixel_adr] = 0;
          for (int block_id = 0; block_id < n_blocks; ++block_id) {

            rho[i * scanner.n_y_pixels + j] +=
                output_rho[block_id * scanner.n_y_pixels + pixel_adr];
          }
        }
      }

#else
      cudaMemcpy(rho.host_ptr,
                 output_rho.device_ptr,
                 image_size,
                 cudaMemcpyDeviceToHost);
#endif
      progress(ib * n_iterations_in_block + it, true);

      // always output first 5 iterations, and at 10, 15, 20, 30, 50, 100
      if (!ib && it < n_iterations_in_block - 1 &&
          (it < 5 || it == 9 || it == 14 || it == 19 || it == 29 || it == 49 ||
           it == 99)) {
        output(it + 1, rho.host_ptr);
      }
    }

    output((ib + 1) * n_iterations_in_block, rho.host_ptr);
  }

  progress(n_iteration_blocks * n_iterations_in_block, false);
}

template <typename F>
void fill_with_sensitivity(F* sensitivity, Scanner<F, short>& scanner) {

  size_t width = scanner.n_z_pixels;
  size_t height = scanner.n_y_pixels;

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      sensitivity[y * width + x] =
          scanner.pixel_sensitivity(Pixel<short>(x, y));
    }
  }
}

// Explicit template instantiation
template void run(Scanner<float, short>& scanner,
                  Response<float>* responses,
                  int n_responses,
                  int n_iteration_blocks,
                  int n_iterations_in_block,
                  util::delegate<void(int, float*)> output,
                  util::delegate<void(int, bool)> progress,
                  int device,
                  int n_blocks,
                  int n_threads_per_block,
                  util::delegate<void(const char*)> device_name);

}  // Reconstruction
}  // GPU
}  // Strip
}  // PET2D
