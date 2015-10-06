#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

#include "util/cuda/debug.h"  // catches all CUDA errors
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
         int n_threads_per_block) {

  cudaSetDevice(device);
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device);

#if NO_EXTRA_OUTPUT
  if (verbose) {
    fprintf(stdout, "# running on: %s\n", prop.name);
  }
#endif

#if __CUDACC__
  dim3 blocks(n_blocks);
  dim3 threads(n_threads_per_block);
#endif

  size_t image_size = scanner.total_n_pixels * sizeof(F);

  const int width = scanner.n_z_pixels;
  const int height = scanner.n_y_pixels;

  cudaChannelFormatDesc desc = cudaCreateChannelDesc<float>();

  F* cpu_sensitivity = new F[scanner.total_n_pixels];

  fill_with_sensitivity(cpu_sensitivity, scanner);

  output(-1, cpu_sensitivity);

  F* gpu_sensitivity;
  size_t pitch_sensitivity;
  cudaMallocPitch(
      &gpu_sensitivity, &pitch_sensitivity, sizeof(F) * width, height);
  cudaMemcpy2D(gpu_sensitivity,
               pitch_sensitivity,
               cpu_sensitivity,
               sizeof(F) * width,
               sizeof(F) * width,
               height,
               cudaMemcpyHostToDevice);
  delete[] cpu_sensitivity;

  cudaBindTexture2D(NULL,
                    &tex_sensitivity,
                    gpu_sensitivity,
                    &desc,
                    width,
                    height,
                    pitch_sensitivity);

  F* cpu_rho = new F[scanner.total_n_pixels];
  for (int i = 0; i < scanner.total_n_pixels; ++i) {
    cpu_rho[i] = 100;
  }

  // this class allocated CUDA pointers and deallocated them in destructor
  GPU::ResponsesSOA<F> gpu_responses(responses, n_responses);

  F* gpu_rho;
  size_t pitch_rho;
  cudaMallocPitch(&gpu_rho, &pitch_rho, sizeof(F) * width, height);
  cudaBindTexture2D(NULL, &tex_rho, gpu_rho, &desc, width, height, pitch_rho);

  F* gpu_output_rho;

#if USE_RHO_PER_WARP
  cudaMalloc((void**)&gpu_output_rho, n_blocks * image_size);
  F* cpu_output_rho;
  cpu_output_rho = new F[n_blocks * scanner.total_n_pixels];
#else
  cudaMalloc((void**)&gpu_output_rho, image_size);
#endif

  for (int ib = 0; ib < n_iteration_blocks; ++ib) {
    for (int it = 0; it < n_iterations_in_block; ++it) {
      progress(ib * n_iterations_in_block + it, false);

#if USE_RHO_PER_WARP
      cudaMemset(gpu_output_rho, 0, n_blocks * image_size);
#else
      cudaMemset(gpu_output_rho, 0, image_size);
#endif
      cudaMemcpy2D(gpu_rho,
                   pitch_rho,
                   cpu_rho,
                   sizeof(F) * width,
                   sizeof(F) * width,
                   height,
                   cudaMemcpyHostToDevice);

#if __CUDACC__
#define reconstruction reconstruction<Kernel><<<blocks, threads>>>
#endif
      reconstruction(scanner,
                     gpu_responses.z_u,
                     gpu_responses.z_d,
                     gpu_responses.dl,
                     n_responses,
                     gpu_output_rho,
                     n_blocks,
                     n_threads_per_block);

      cudaThreadSynchronize();

#if USE_RHO_PER_WARP
      cudaMemcpy(cpu_output_rho,
                 gpu_output_rho,
                 n_blocks * image_size,
                 cudaMemcpyDeviceToHost);

      for (int i = 0; i < scanner.n_y_pixels; ++i) {
        for (int j = 0; j < scanner.n_z_pixels; ++j) {
          int pixel_adr = i * scanner.n_y_pixels + j;
          cpu_rho[pixel_adr] = 0;
          for (int block_id = 0; block_id < n_blocks; ++block_id) {

            cpu_rho[i * scanner.n_y_pixels + j] +=
                cpu_output_rho[block_id * scanner.n_y_pixels + pixel_adr];
          }
        }
      }

#else
      cudaMemcpy(cpu_rho, gpu_output_rho, image_size, cudaMemcpyDeviceToHost);
#endif
      progress(ib * n_iterations_in_block + it, true);

      // always output first 5 iterations, and at 10, 15, 20, 30, 50, 100
      if (!ib && it < n_iterations_in_block - 1 &&
          (it < 5 || it == 9 || it == 14 || it == 19 || it == 29 || it == 49 ||
           it == 99)) {
        output(it + 1, cpu_rho);
      }
    }

    output((ib + 1) * n_iterations_in_block, cpu_rho);
  }

  progress(n_iteration_blocks * n_iterations_in_block, false);

  cudaUnbindTexture(&tex_sensitivity);
  cudaFree(gpu_sensitivity);
  cudaUnbindTexture(&tex_rho);
  cudaFree(gpu_rho);
  cudaFree(gpu_output_rho);
  delete[] cpu_rho;
#if USE_RHO_PER_WARP
  delete[] cpu_output_rho;
#endif
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
                  int n_threads_per_block);

}  // Reconstruction
}  // GPU
}  // Strip
}  // PET2D
