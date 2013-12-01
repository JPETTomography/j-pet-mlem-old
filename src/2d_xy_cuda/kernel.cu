// if we don't include that Qt Creator will show many errors
#include <cuda_runtime.h>
#include <sys/time.h>
#include <stdio.h>
#include "config.h"
#include "data_structures.h"
#include "prng.cuh"
#include "geometry_methods.cuh"
#include "gpu_detector_geometry_test.cuh"
#include "gpu_detector_hit_test.cuh"
#include "gpu_phantom_generation.cuh"

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

void mem_clean_lors(Matrix_Element* cpu_matrix, int number_of_blocks) {

  for (int i = 0; i < number_of_blocks; ++i) {
    for (int j = 0; j < LORS; ++j) {

      cpu_matrix[i].lor[j] = 0.f;
    }
  }
}

void phantom_kernel(int number_of_threads_per_block,
                    int number_of_blocks,
                    int n_emissions,
                    int pixels_in_row,
                    float radius,
                    float h_detector,
                    float w_detector,
                    float pixel_size) {

  dim3 blocks(number_of_blocks);
  dim3 threads(number_of_threads_per_block);

  unsigned int* cpu_prng_seed;

  cudaSetDevice(0);

  cpu_prng_seed =
      (unsigned int*)malloc(number_of_blocks * number_of_threads_per_block * 4 *
                            sizeof(unsigned int));

  for (int i = 0; i < 4 * number_of_blocks * number_of_threads_per_block; ++i) {

    cpu_prng_seed[i] = 53445 + i;
  }

  int triangular_matrix_size =
      ((pixels_in_row / 2) * ((pixels_in_row / 2) + 1) / 2);

  Matrix_Element* cpu_matrix =
      (Matrix_Element*)malloc(number_of_blocks * sizeof(Matrix_Element));

  // unsigned int matrix_size = triangular_matrix_size * number_of_blocks;

  unsigned int* gpu_prng_seed;
  Matrix_Element* gpu_matrix_element;

  cuda(Malloc,
       (void**)&gpu_prng_seed,
       number_of_blocks * number_of_threads_per_block * 4 *
           sizeof(unsigned int));
  cuda(Malloc,
       (void**)&gpu_matrix_element,
       number_of_blocks * sizeof(Matrix_Element));

  cuda(
      Memcpy,
      gpu_prng_seed,
      cpu_prng_seed,
      number_of_blocks * number_of_threads_per_block * 4 * sizeof(unsigned int),
      cudaMemcpyHostToDevice);

  printf("GPU kernel start\n");
  printf(
      "Number of Detectors %d Numer of LORS: %d\n", NUMBER_OF_DETECTORS, LORS);

  double timer = getwtime();

  //  for (int j = pixels_in_row / 2 - 1; j >= 0; --j) {
  //    for (int i = 0; i <= j; ++i) {

  for (int j = 0; j < 1; ++j) {
    for (int i = 0; i < 1; ++i) {

      mem_clean_lors(cpu_matrix, number_of_blocks);

      cuda(Memcpy,
           gpu_matrix_element,
           cpu_matrix,
           number_of_blocks * sizeof(Matrix_Element),
           cudaMemcpyHostToDevice);

      printf("Pixel(%d,%d) n_emissions: %d \n", i, j, n_emissions);

      gpu_phantom_generation << <blocks, threads>>>
          (i,
           j,
           n_emissions,
           gpu_prng_seed,
           gpu_matrix_element,
           number_of_threads_per_block,
           pixels_in_row,
           radius,
           h_detector,
           w_detector,
           pixel_size);

      cudaThreadSynchronize();

      cuda(Memcpy,
           cpu_matrix,
           gpu_matrix_element,
           number_of_blocks * sizeof(Matrix_Element),
           cudaMemcpyDeviceToHost);

      for (int i = 0; i < LORS; i++) {
        float temp = 0.f;
        for (int j = 1; j < 2; ++j) {

          temp += cpu_matrix[j].lor[i];
        }

        if (temp > 0.0f) {
          printf("%f\n",
                 temp / number_of_blocks / number_of_threads_per_block /
                     n_emissions);
        }
      }
    }
  }
  double time = 0.0f;

  time = getwtime() - time;

  printf("time[s]: %f\n ", time);
  printf("time per pixel: %f\n", time / triangular_matrix_size);

  cuda(Free, gpu_prng_seed);
  cuda(Free, gpu_matrix_element);
}

void gpu_detector_geometry_kernel_test(float radius,
                                       float h_detector,
                                       float w_detector,
                                       float pixel_size,
                                       Detector_Ring& cpu_output) {

  dim3 blocks(1);
  dim3 threads(NUMBER_OF_DETECTORS);

  cudaSetDevice(0);

  Detector_Ring* cpu_detectors = (Detector_Ring*)malloc(sizeof(Detector_Ring));

  Detector_Ring* gpu_detectors;

  cuda(Malloc, (void**)&gpu_detectors, sizeof(Detector_Ring));

  printf("Execute gpu_kernel_test for detectors geometry\n");

  cuda(Memcpy,
       gpu_detectors,
       cpu_detectors,
       sizeof(Detector_Ring),
       cudaMemcpyHostToDevice);

  gpu_detector_geometry_test << <blocks, threads>>>
      (radius, h_detector, w_detector, pixel_size, gpu_detectors);

  cudaThreadSynchronize();

  cuda(Memcpy,
       cpu_detectors,
       gpu_detectors,
       sizeof(Detector_Ring),
       cudaMemcpyDeviceToHost);

  for (int i = 0; i < NUMBER_OF_DETECTORS; ++i) {
    for (int j = 0; j < 4; j++) {

      cpu_output.detector_list[i].points[j].x =
          cpu_detectors->detector_list[i].points[j].x;
      cpu_output.detector_list[i].points[j].y =
          cpu_detectors->detector_list[i].points[j].y;
    }
  }

  cuda(Free, gpu_detectors);
}

void gpu_detector_hits_kernel_test(float crx,
                                   float cry,
                                   float cangle,
                                   float radius,
                                   float h_detector,
                                   float w_detector) {

  cudaSetDevice(0);

  gpu_detector_hit_test << <1, NUMBER_OF_DETECTORS>>>
      (crx, cry, cangle, radius, h_detector, w_detector);

  cudaThreadSynchronize();
}
