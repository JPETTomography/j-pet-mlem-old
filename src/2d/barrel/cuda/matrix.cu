#include <cuda_runtime.h>
#include <stdio.h>

#include "util/cuda/debug.h"  // catches all CUDA errors

#include "config.h"
#include "prng.cuh"
#include "geometry.h"
#include "monte_carlo.cuh"

#include "2d/geometry/pixel.h"

namespace PET2D {
namespace Barrel {
namespace GPU {

bool run_matrix(Pixel<> pixel,
                int n_tof_positions,
                int number_of_threads_per_block,
                int number_of_blocks,
                int n_emissions,
                float radius,
                float h_detector,
                float w_detector,
                float pixel_size,
                unsigned int* cpu_prng_seed,
                MatrixElement* cpu_matrix,
                MatrixElement* gpu_output) {

  dim3 blocks(number_of_blocks);
  dim3 threads(number_of_threads_per_block);

  cudaSetDevice(0);

  unsigned int* gpu_prng_seed;
  MatrixElement* gpu_MatrixElement;

#if WARP_DIVERGENCE_TEST
  bool* warp_divergence_buffer;

  const int warp_size = 32;

  cudaMalloc((void**)&warp_divergence_buffer,
             warp_size * n_emissions * sizeof(bool));

  cudaMemset(warp_divergence_buffer, 0, warp_size * n_emissions * sizeof(bool));

#else
  bool* warp_divergence_buffer;
  cudaMalloc((void**)&warp_divergence_buffer, 0 * sizeof(bool));

#endif

  cudaMalloc((void**)&gpu_prng_seed,
             number_of_blocks * number_of_threads_per_block * 4 *
                 sizeof(unsigned int));
  cudaMalloc((void**)&gpu_MatrixElement,
             n_tof_positions * number_of_blocks * sizeof(MatrixElement));

  cudaMemcpy(
      gpu_prng_seed,
      cpu_prng_seed,
      number_of_blocks * number_of_threads_per_block * 4 * sizeof(unsigned int),
      cudaMemcpyHostToDevice);

  float fov_radius = radius / M_SQRT2;

  int i = pixel.x;
  int j = pixel.y;

  cudaMemset(gpu_MatrixElement,
             0,
             n_tof_positions * number_of_blocks * sizeof(MatrixElement));

  long total_emissions =
      (long)n_emissions * number_of_blocks * number_of_threads_per_block;

  printf(
      "Pixel(%d,%d) n_emissions: %d %ld\n", i, j, n_emissions, total_emissions);

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start);

  if ((i * i + j * j) * pixel_size * pixel_size < fov_radius * fov_radius) {

#if __CUDACC__
#define monte_carlo_kernel monte_carlo_kernel << <blocks, threads>>>
#endif
    monte_carlo_kernel(i,
                       j,
                       n_emissions,
                       n_tof_positions,
                       gpu_prng_seed,
                       gpu_MatrixElement,
                       radius,
                       h_detector,
                       w_detector,
                       pixel_size,
                       warp_divergence_buffer);

    cudaThreadSynchronize();

    if (cudaGetLastError() != cudaSuccess) {
      return false;
    }
  }

  cudaEventRecord(stop);

  cudaEventSynchronize(stop);
  float milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);

  printf("Direct kernel time without memcpy %f ms\n", milliseconds);

#if WARP_DIVERGENCE_TEST

  bool* cpu_warp_divergence_buffer = new bool[warp_size * n_emissions];

  cuda(Memcpy,
       cpu_warp_divergence_buffer,
       warp_divergence_buffer,
       warp_size * n_emissions * sizeof(bool),
       cudaMemcpyDeviceToHost);

  std::ofstream output;
  output.open("warp_info");

  for (int i = 0; i < warp_size * n_emissions; ++i) {

    output << int(cpu_warp_divergence_buffer[i]);
    if (i % warp_size == 0 && i != 0) {
      output << std::endl;
    }
  }

  delete cpu_warp_divergence_buffer;

#endif

  cudaMemcpy(cpu_matrix,
             gpu_MatrixElement,
             n_tof_positions * number_of_blocks * sizeof(MatrixElement),
             cudaMemcpyDeviceToHost);

#if NO_TOF > 0
  for (int lor_i = 0; lor_i < LORS; ++lor_i) {
    float temp = 0;
    for (int block_i = 0; block_i < number_of_blocks; ++block_i) {

      temp += cpu_matrix[block_i].hit[lor_i];
    }

    if (temp > 0) {
      gpu_output[0].hit[lookup_table_lors[lor_i].index()] = temp;
    }
  }
#else
  for (int tof_i = 0; tof_i < n_tof_positions; ++tof_i) {
    for (int lor_i = 0; lor_i < LORS; ++lor_i) {
      float temp_hits = 0;
      for (int block_i = 0; block_i < number_of_blocks; ++block_i) {
        temp_hits += cpu_matrix[tof_i + (block_i * n_tof_positions)].hit[lor_i];
      }
      if (temp_hits > 0) {
        gpu_output[tof_i].hit[lor_i] = temp_hits;
      }
    }
  }

#endif
  cudaFree(gpu_prng_seed);
  cudaFree(gpu_MatrixElement);

  return 0;
}

}  // GPU
}  // Barrel
}  // PET2D
