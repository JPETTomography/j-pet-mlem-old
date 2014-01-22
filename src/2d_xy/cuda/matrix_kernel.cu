#include <cuda_runtime.h>

#include <sys/time.h>

#include "config.h"
#include "prng.cuh"
#include "geometry.h"
#include "geometry_methods.cuh"
#include "detector_geometry_test.cuh"
#include "detector_hit_test.cuh"
#include "monte_carlo.cuh"

// FIXME: remove me
#include "geometry/pixel.h"

using namespace gpu;

static cudaError err;

#define cuda(f, ...)                                        \
  if ((err = cuda##f(__VA_ARGS__)) != cudaSuccess) {        \
    fprintf(stderr, #f "() %s\n", cudaGetErrorString(err)); \
    exit(-1);                                               \
  }
#define cudathread_per_blockoSync(...) cuda(__VA_ARGS__)

double getwtime() {
  struct timeval tv;
  static time_t sec = 0;
  gettimeofday(&tv, NULL);
  if (!sec)
    sec = tv.tv_sec;
  return (double)(tv.tv_sec - sec) + (double)tv.tv_usec / 1e6;
}

void mem_clean_lors(MatrixElement* cpu_matrix, int number_of_blocks) {

  for (int i = 0; i < number_of_blocks; ++i) {
    for (int j = 0; j < LORS; ++j) {

      cpu_matrix[i].hit[j] = 0.f;
    }
  }
}

bool run_monte_carlo_kernel(int pixel_i,
                            int number_of_threads_per_block,
                            int number_of_blocks,
                            int n_emissions,
                            int n_detectors,
                            int pixels_in_row,
                            int triangular_matrix_size,
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

  cudaSetDevice(1);

  unsigned int* gpu_prng_seed;
  MatrixElement* gpu_MatrixElement;

  cuda(Malloc,
       (void**)&gpu_prng_seed,
       number_of_blocks * number_of_threads_per_block * 4 *
           sizeof(unsigned int));
  cuda(Malloc,
       (void**)&gpu_MatrixElement,
       number_of_blocks * sizeof(MatrixElement));

  cuda(
      Memcpy,
      gpu_prng_seed,
      cpu_prng_seed,
      number_of_blocks * number_of_threads_per_block * 4 * sizeof(unsigned int),
      cudaMemcpyHostToDevice);

  float fov_radius = radius / M_SQRT2;

  unsigned int number_of_lors = (n_detectors * (n_detectors + 1)) / 2;

  Pixel<> pixel = lookup_table_pixel[pixel_i];

  int i = pixel.x;
  int j = pixel.y;

  mem_clean_lors(cpu_matrix, number_of_blocks);

  cuda(Memcpy,
       gpu_MatrixElement,
       cpu_matrix,
       number_of_blocks * sizeof(MatrixElement),
       cudaMemcpyHostToDevice);

  long total_emissions =
      (long)n_emissions * number_of_blocks * number_of_threads_per_block;

  printf(
      "Pixel(%d,%d) n_emissions: %d %ld\n", i, j, n_emissions, total_emissions);

  if ((i * i + j * j) * pixel_size * pixel_size < fov_radius * fov_radius) {

    monte_carlo_kernel << <blocks, threads>>> (i,
                                               j,
                                               n_emissions,
                                               n_detectors,
                                               gpu_prng_seed,
                                               gpu_MatrixElement,
                                               number_of_threads_per_block,
                                               pixels_in_row,
                                               radius,
                                               h_detector,
                                               w_detector,
                                               pixel_size);

    cudaThreadSynchronize();

    if (cudaGetLastError() != cudaSuccess) {
      return false;
    }
  }

  cuda(Memcpy,
       cpu_matrix,
       gpu_MatrixElement,
       number_of_blocks * sizeof(MatrixElement),
       cudaMemcpyDeviceToHost);

  for (int i = 0; i < LORS; i++) {
    float temp = 0.f;
    for (int j = 0; j < number_of_blocks; ++j) {

      temp += cpu_matrix[j].hit[i];
    }

    if (temp > 0.0f) {
      gpu_output[0].hit[lookup_table_lors[i].index()] = temp;
    }
  }

  cuda(Free, gpu_prng_seed);
  cuda(Free, gpu_MatrixElement);

  return 0;
}
