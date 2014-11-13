// FIXME: stupid workaround for GCC or NVCC bug causing errors if this not gets
// included first
#include <type_traits>

#include <cuda_runtime.h>
#include <cstdio>
#include <ctime>

#include "util/cuda/debug.h"  // catches all CUDA errors

#include "phantom.h"

namespace PET2D {
namespace Barrel {
namespace GPU {

void run_gpu_phantom(int n_threads_per_block,
                     int n_blocks,
                     int n_emissions,
                     int n_pixels_in_row,
                     float s_size,
                     Pixel* lookup_table_pixel,
                     int* pixel_hits) {
  cudaSetDevice(1);
#if PHANTOM_IS_BROKEN_FOR_NOW
  dim3 blocks(n_blocks);
  dim3 threads(n_threads_per_block);

  unsigned int* cpu_prng_seed;

  cpu_prng_seed = new unsigned int[n_blocks * n_threads_per_block * 4];

  for (int i = 0; i < 4 * n_blocks * n_threads_per_block; ++i) {

    cpu_prng_seed[i] = 53445 + i;
  }

  int triangular_matrix_size =
      ((n_pixels_in_row / 2) * ((n_pixels_in_row / 2) + 1) / 2);

  for (int i = 0; i < triangular_matrix_size; ++i) {

    for (int lor = 0; lor < LORS; ++lor) {

      gpu_output[i].hit[lor] = 0;
    }
  }

  MatrixElement* cpu_matrix = new MatrixElement[n_blocks];

  unsigned int* gpu_prng_seed;
  MatrixElement* gpu_MatrixElement;

  cudaMalloc((void**)&gpu_prng_seed,
             n_blocks * n_threads_per_block * 4 * sizeof(unsigned int));
  cudaMalloc((void**)&gpu_MatrixElement, n_blocks * sizeof(MatrixElement));

  cudaMemcpy(gpu_prng_seed,
             cpu_prng_seed,
             n_blocks * n_threads_per_block * 4 * sizeof(unsigned int),
             cudaMemcpyHostToDevice);

  printf("GPU kernel start\n");
  printf("DETECTORS %d LORS: %d\n", NUMBER_OF_DETECTORS, LORS);

  for (int p = 0; p < triangular_matrix_size; ++p) {

    Pixel pixel = lookup_table_pixel[p];

    int i = pixel.x;
    int j = pixel.y;
#if BROKEN
    mem_clean_lors(cpu_matrix, number_of_blocks);
#endif
    cudaMemcpy(gpu_MatrixElement,
               cpu_matrix,
               n_blocks * sizeof(MatrixElement),
               cudaMemcpyHostToDevice);

    long total_emissions = (long)n_emissions * n_blocks * n_threads_per_block;

    printf("Pixel(%d,%d) n_emissions: %d %ld\n",
           i,
           j,
           n_emissions,
           total_emissions);

    float fov_radius = radius / M_SQRT2;
    if ((i * i + j * j) * s_size * s_size < fov_radius * fov_radius) {
#if __CUDACC__
#define gpu_phantom_generation gpu_phantom_generation << <blocks, threads>>>
#endif
#if WHERE_IS_PHANTOM
      gpu_phantom_generation(i,
                             j,
                             n_emissions,
                             gpu_prng_seed,
                             gpu_MatrixElement,
                             number_of_threads_per_block,
                             pixels_in_row,
                             radius,
                             h_detector,
                             w_detector,
                             pixel_size);
#endif
      cudaThreadSynchronize();
    }
    cudaMemcpy(cpu_matrix,
               gpu_MatrixElement,
               n_blocks * sizeof(MatrixElement),
               cudaMemcpyDeviceToHost);

    if (p == 0) {
      for (int i = 0; i < LORS; i++) {
        float temp = 0.f;
        for (int j = 0; j < n_blocks; ++j) {

          temp += cpu_matrix[j].hit[i];
        }
#if BROKEN
        if (temp > 0.0f) {
          GPU::LOR lor(lookup_table_lors[i].lor_a, lookup_table_lors[i].lor_b);
          gpu_output[p].hit[lor.index()] = temp;
        }
#endif
      }
    }
  }

  cudaFree(gpu_prng_seed);
  cudaFree(gpu_MatrixElement);
#endif
}

}  // GPU
}  // Barrel
}  // PET2D
