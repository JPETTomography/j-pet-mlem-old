#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "util/cuda/debug.h"  // catches all CUDA errors
#include "../event.h"
#include "../kernel.h"
#include "gpu_events_soa.h"
#include "config.h"

texture<float, 2, cudaReadModeElementType> tex_sensitivity;
texture<float, 2, cudaReadModeElementType> tex_rho;

#if THREAD_GRANULARITY
#include "reconstruction_thread_granularity.cuh"
#elif WARP_GRANULARITY
#include "reconstruction_warp_granularity.cuh"
#else
#include "reconstruction_simple.cuh"
#endif

template <typename F>
void fill_with_sensitivity(F* sensitivity,
                           F* inv_sensitivity,
                           StripDetector<F>& detector);

template <typename F>
void run_gpu_reconstruction(StripDetector<F>& detector,
                            Event<F>* events,
                            int n_events,
                            int n_iteration_blocks,
                            int n_iterations_in_block,
                            void (*output_callback)(StripDetector<F>& detector,
                                                    int iteration,
                                                    F* image,
                                                    void* context),
                            void (*progress_callback)(int iteration,
                                                      void* context),
                            void* context,
                            int device,
                            int n_blocks,
                            int n_threads_per_block,
                            bool verbose) {

  cudaSetDevice(device);
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device);

  if (verbose) {
    fprintf(stdout, "Running on: %s\n", prop.name);
  }

#if __CUDACC__
  dim3 blocks(n_blocks);
  dim3 threads(n_threads_per_block);
#endif

  size_t image_size = detector.total_n_pixels * sizeof(F);

  const int width = detector.n_z_pixels;
  const int height = detector.n_y_pixels;

  cudaChannelFormatDesc desc = cudaCreateChannelDesc<float>();

  F* cpu_sensitivity = new F[detector.total_n_pixels];

  fill_with_sensitivity(cpu_sensitivity, detector);

  output_callback(detector, -1, cpu_sensitivity, context);

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

  F* cpu_rho = new F[detector.total_n_pixels];
  for (int i = 0; i < detector.total_n_pixels; ++i) {
    cpu_rho[i] = 100;
  }

  // this class allocated CUDA pointers and deallocated them in destructor
  GPUEventsSOA<F> gpu_events(events, n_events);

  F* gpu_rho;
  size_t pitch_rho;
  cudaMallocPitch(&gpu_rho, &pitch_rho, sizeof(F) * width, height);
  cudaBindTexture2D(NULL, &tex_rho, gpu_rho, &desc, width, height, pitch_rho);

  F* gpu_output_rho;

#if USE_RHO_PER_WARP
  cudaMalloc((void**)&gpu_output_rho, n_blocks * image_size);
  F* cpu_output_rho;
  cpu_output_rho = new F[n_blocks * detector.total_n_pixels];
#else
  cudaMalloc((void**)&gpu_output_rho, image_size);
#endif

  for (int ib = 0; ib < n_iteration_blocks; ++ib) {
    for (int it = 0; it < n_iterations_in_block; ++it) {

      cudaEvent_t start, stop, start_mem_time, stop_mem_time;
      float time;
      float time_all;
      if (verbose) {
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        cudaEventCreate(&start_mem_time);
        cudaEventCreate(&stop_mem_time);
      } else {
        progress_callback(ib * n_iterations_in_block + it, context);
      }

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

      if (verbose) {
        cudaEventRecord(start);
        cudaEventRecord(start_mem_time);
      }
#ifdef __METRIC__
      F* gpu_metric_memory;
      const int metric_size = n_blocks * n_threads_per_block;
      cudaMalloc((void**)&gpu_metric_memory, metric_size);
      cudaMemset(gpu_metric_memory, 0, metric_size);
#endif

#if __CUDACC__
#define reconstruction reconstruction<Kernel> << <blocks, threads>>>
#endif
#ifdef __METRIC__
      reconstruction(gpu_metric_memory,
                     detector,
                     gpu_events.z_u,
                     gpu_events.z_d,
                     gpu_events.dl,
                     n_events,
                     gpu_output_rho,
                     n_blocks,
                     n_threads_per_block);
#else
      reconstruction(detector,
                     gpu_events.z_u,
                     gpu_events.z_d,
                     gpu_events.dl,
                     n_events,
                     gpu_output_rho,
                     n_blocks,
                     n_threads_per_block);
#endif

      cudaThreadSynchronize();

      if (verbose) {
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&time, start, stop);
      }

#ifdef __METRIC__
  cudaFree(gpu_metric_memory);
#endif

#if USE_RHO_PER_WARP
      cudaMemcpy(cpu_output_rho,
                 gpu_output_rho,
                 n_blocks * image_size,
                 cudaMemcpyDeviceToHost);

      for (int i = 0; i < detector.n_y_pixels; ++i) {
        for (int j = 0; j < detector.n_z_pixels; ++j) {
          int pixel_adr = i * detector.n_y_pixels + j;
          cpu_rho[pixel_adr] = 0;
          for (int block_id = 0; block_id < n_blocks; ++block_id) {

            cpu_rho[i * detector.n_y_pixels + j] +=
                cpu_output_rho[block_id * detector.n_y_pixels + pixel_adr];
          }
        }
      }

#else
      cudaMemcpy(cpu_rho, gpu_output_rho, image_size, cudaMemcpyDeviceToHost);
#endif

      if (verbose) {
        cudaEventRecord(stop_mem_time);
        cudaEventSynchronize(stop_mem_time);
        cudaEventElapsedTime(&time_all, start_mem_time, stop_mem_time);
        printf(
            "[%02d] kernel       : %f ms\n"
            "     kernel + mem : %f ms\n",
            ib * n_iterations_in_block + it,
            time,
            time_all);
      }
    }

    output_callback(detector, ib * n_iterations_in_block, cpu_rho, context);
  }

  if (!verbose) {
    progress_callback(n_iteration_blocks * n_iterations_in_block, context);
  }

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
void fill_with_sensitivity(F* sensitivity, StripDetector<F>& detector) {

  size_t width = detector.n_z_pixels;
  size_t height = detector.n_y_pixels;

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      sensitivity[y * width + x] = detector.pixel_sensitivity(Pixel<>(x, y));
    }
  }
}

// Explicit template instantiation

template void run_gpu_reconstruction<float>(
    StripDetector<float>& detector,
    Event<float>* events,
    int n_events,
    int n_iteration_blocks,
    int n_iterations_in_block,
    void (*output_callback)(StripDetector<float>& detector,
                            int iteration,
                            float* image,
                            void* context),
    void (*progress_callback)(int iteration, void* context),
    void* context,
    int device,
    int n_blocks,
    int n_threads_per_block,
    bool verbose);
