#include <cuda_runtime.h>
#include <stdio.h>

#include "config.h"
#include "prng.cuh"
#include "geometry.h"
#include "geometry_methods.cuh"
#include "monte_carlo.cuh"

#include "geometry/pixel.h"

using namespace gpu;

static cudaError err;

#define cuda(f, ...)                                        \
  if ((err = cuda##f(__VA_ARGS__)) != cudaSuccess) {        \
    fprintf(stderr, #f "() %s\n", cudaGetErrorString(err)); \
    exit(-1);                                               \
  }
#define cudathread_per_blockoSync(...) cuda(__VA_ARGS__)

bool run_monte_carlo_kernel(int pixel_i,
                            int n_tof_positions,
                            int number_of_threads_per_block,
                            int number_of_blocks,
                            int n_emissions,
                            float radius,
                            float h_detector,
                            float w_detector,
                            float pixel_size,
                            gpu::LOR* lookup_table_lors,
                            Pixel<>* lookup_table_pixel,
                            unsigned int* cpu_prng_seed,
                            MatrixElement* cpu_matrix,
                            MatrixElement* gpu_output) {

  dim3 blocks(number_of_blocks);
  dim3 threads(number_of_threads_per_block);

  cudaSetDevice(0);

  unsigned int* gpu_prng_seed;
  MatrixElement* gpu_MatrixElement;

#if WARP_DIVERGENCE_TEST
  bool* warp_divergence_buffor;

  const int warp_size = 32;

  cuda(Malloc,
       (void**)&warp_divergence_buffor,
       warp_size * n_emissions * sizeof(bool));

  cuda(Memset,
       warp_divergence_buffor,
       0,
       warp_size * n_emissions * sizeof(bool));

#else
  bool* warp_divergence_buffor;
  cuda(Malloc, (void**)&warp_divergence_buffor, 0 * sizeof(bool));

#endif

  cuda(Malloc,
       (void**)&gpu_prng_seed,
       number_of_blocks * number_of_threads_per_block * 4 *
           sizeof(unsigned int));
  cuda(Malloc,
       (void**)&gpu_MatrixElement,
       n_tof_positions * number_of_blocks * sizeof(MatrixElement));

  cuda(
      Memcpy,
      gpu_prng_seed,
      cpu_prng_seed,
      number_of_blocks * number_of_threads_per_block * 4 * sizeof(unsigned int),
      cudaMemcpyHostToDevice);

  float fov_radius = radius / M_SQRT2;

  Pixel<> pixel = lookup_table_pixel[pixel_i];

  int i = pixel.x;
  int j = pixel.y;

  cuda(Memset,
       gpu_MatrixElement,
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

    monte_carlo_kernel << <blocks, threads>>> (i,
                                               j,
                                               n_emissions,
                                               n_tof_positions,
                                               gpu_prng_seed,
                                               gpu_MatrixElement,
                                               radius,
                                               h_detector,
                                               w_detector,
                                               pixel_size,
                                               warp_divergence_buffor);

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

  bool* cpu_warp_divergence_buffor = new bool[warp_size * n_emissions];

  cuda(Memcpy,
       cpu_warp_divergence_buffor,
       warp_divergence_buffor,
       warp_size * n_emissions * sizeof(bool),
       cudaMemcpyDeviceToHost);

  std::ofstream output;
  output.open("warp_info");

  for (int i = 0; i < warp_size * n_emissions; ++i) {

    output << int(cpu_warp_divergence_buffor[i]);
    if (i % warp_size == 0 && i != 0) {
      output << std::endl;
    }
  }

  delete cpu_warp_divergence_buffor;

#endif

  cuda(Memcpy,
       cpu_matrix,
       gpu_MatrixElement,
       n_tof_positions * number_of_blocks * sizeof(MatrixElement),
       cudaMemcpyDeviceToHost);

#if NO_TOF > 0
  for (int lor_i = 0; lor_i < LORS; ++lor_i) {
    float temp = 0.f;
    for (int block_i = 0; block_i < number_of_blocks; ++block_i) {

      temp += cpu_matrix[block_i].hit[lor_i];
    }

    if (temp > 0.0f) {
      gpu_output[0].hit[lookup_table_lors[lor_i].index()] = temp;
    }
  }
#else
  for (int tof_i = 0; tof_i < n_tof_positions; ++tof_i) {
    for (int lor_i = 0; lor_i < LORS; ++lor_i) {
      float temp_hits = 0.f;
      for (int block_i = 0; block_i < number_of_blocks; ++block_i) {

        temp_hits += cpu_matrix[tof_i + (block_i * n_tof_positions)].hit[lor_i];
      }

      if (temp_hits > 0.0f) {

        gpu_output[tof_i].hit[lookup_table_lors[lor_i].index()] = temp_hits;
      }
    }
  }

#endif
  cuda(Free, gpu_prng_seed);
  cuda(Free, gpu_MatrixElement);

  return 0;
}
