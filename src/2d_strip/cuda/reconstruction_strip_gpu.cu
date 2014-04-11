#include <cuda_runtime.h>
#include <iostream>
#include <fstream>
#include "reconstruction.cuh"

static cudaError err;

#define cuda_kernel_config(blocks, threads)                      \
  {                                                              \
    printf("Cuda kernel config\n");                              \
    printf("Number of  blocks per kernel: %d\n", blocks);        \
    printf("Number of threads|block per kernel: %d\n", threads); \
  }

#define cuda(f, ...)                                        \
  if ((err = cuda##f(__VA_ARGS__)) != cudaSuccess) {        \
    fprintf(stderr, #f "() %s\n", cudaGetErrorString(err)); \
    exit(-1);                                               \
  }
#define cudathread_per_blockoSync(...) cuda(__VA_ARGS__)

void gpu_reconstruction_strip_2d(gpu_config::GPU_parameters cfg,
                                 event<float>* event_list,
                                 int iteration_chunk) {

  dim3 blocks(cfg.number_of_blocks);
  dim3 threads(cfg.number_of_threads_per_block);

  cuda_kernel_config(cfg.number_of_blocks, cfg.number_of_threads_per_block);

  cudaSetDevice(0);

  int* cpu_example_data;

  cpu_example_data = (int*)malloc(1000 * sizeof(int));

  int* gpu_example_data;

  cuda(Malloc, (void**)&gpu_example_data, 1000 * sizeof(int));

  cuda(Memcpy,
       gpu_example_data,
       cpu_example_data,
       1000 * sizeof(int),
       cudaMemcpyHostToDevice);

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start);

  reconstruction_2d_strip_cuda<float> << <blocks, threads>>>
      (cfg, event_list, iteration_chunk);

  cudaThreadSynchronize();

  //  if (cudaGetLastError() != cudaSuccess) {
  //    return false;
  //  }

  cudaEventRecord(stop);

  cudaEventSynchronize(stop);
  float milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);

  printf("Direct kernel time without memcpy %f ms\n", milliseconds);

  cuda(Free, gpu_example_data);
  free(cpu_example_data);
}
