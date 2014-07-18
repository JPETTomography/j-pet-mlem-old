#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

#include "util/cuda/debug.h"  // catches all CUDA errors
#include "../event.h"

#include "config.h"

#if USE_TEXTURE
texture<float, 2, cudaReadModeElementType> tex_sensitivity;
#endif

#if EVENT_GRANULARITY
#include "reconstruction_event_granularity.cuh"
#elif WARP_GRANULARITY
#include "reconstruction_warp_granularity.cuh"
#else
#include "reconstruction_simple.cuh"
#endif

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
                            int n_threads_per_block) {

  cudaSetDevice(device);

  dim3 blocks(n_blocks);
  dim3 threads(n_threads_per_block);

  size_t image_size = detector.total_n_pixels * sizeof(F);
  size_t output_size = image_size * n_blocks;
  size_t events_size = n_events * sizeof(F);

  F* cpu_sensitivity = (F*)malloc(image_size);
  F* cpu_output = (F*)malloc(output_size);
  F* cpu_rho = (F*)malloc(image_size);

  for (int y = 0; y < detector.n_y_pixels; ++y) {
    for (int x = 0; x < detector.n_z_pixels; ++x) {
      Point<F> point = detector.pixel_center(y, x);
      cpu_sensitivity[y * detector.n_z_pixels + x] =
          detector.sensitivity(point);
    }
  }

  output_callback(detector, -1, cpu_sensitivity, context);

  for (int i = 0; i < detector.total_n_pixels; ++i) {
    cpu_rho[i] = 100;
  }

  F* cpu_events_z_u = (F*)malloc(events_size);
  F* cpu_events_z_d = (F*)malloc(events_size);
  F* cpu_events_dl = (F*)malloc(events_size);

  for (int i = 0; i < n_events; ++i) {
    cpu_events_z_u[i] = events[i].z_u;
    cpu_events_z_d[i] = events[i].z_d;
    cpu_events_dl[i] = events[i].dl;
  }

  F* gpu_sensitivity;
  F* gpu_output;
  F* gpu_rho;

  size_t pitch;
  cudaMallocPitch(&gpu_sensitivity,
                  &pitch,
                  sizeof(F) * detector.n_z_pixels,
                  detector.n_y_pixels);

  cudaMemcpy2D(gpu_sensitivity,
               pitch,
               cpu_sensitivity,
               sizeof(F) * detector.n_z_pixels,
               sizeof(F) * detector.n_z_pixels,
               detector.n_y_pixels,
               cudaMemcpyHostToDevice);

  free(cpu_sensitivity);

#if USE_TEXTURE_OBJECT
  // create texture description
  cudaResourceDesc resDesc;
  memset(&resDesc, 0, sizeof(resDesc));
  resDesc.resType = cudaResourceTypePitch2D;
  resDesc.res.pitch2D.devPtr = gpu_sensitivity;
  resDesc.res.pitch2D.pitchInBytes = pitch;
  resDesc.res.pitch2D.width = detector.n_y_pixels;
  resDesc.res.pitch2D.height = detector.n_z_pixels;
  // resDesc.res.pitch2D.desc = cudaCreateChannelDesc<F>();
  resDesc.res.pitch2D.desc.f = cudaChannelFormatKindFloat;
  resDesc.res.pitch2D.desc.x = 32;  // 32 bits per channel for float texture
  cudaTextureDesc texDesc;
  memset(&texDesc, 0, sizeof(texDesc));
  texDesc.readMode = cudaReadModeElementType;
  // create texture object
  cudaTextureObject_t tex_sensitivity;
  cudaCreateTextureObject(&tex_sensitivity, &resDesc, &texDesc, NULL);
#elif USE_TEXTURE
  cudaChannelFormatDesc desc_sensitivity = cudaCreateChannelDesc<float>();
  cudaBindTexture2D(NULL,
                    &tex_sensitivity,
                    gpu_sensitivity,
                    &desc_sensitivity,
                    detector.n_z_pixels,
                    detector.n_y_pixels,
                    pitch);
#endif

  cudaMalloc((void**)&gpu_output, output_size);
  cudaMalloc((void**)&gpu_rho, image_size);

  F* gpu_events_z_u;
  F* gpu_events_z_d;
  F* gpu_events_dl;

  cudaMalloc((void**)&gpu_events_z_u, events_size);
  cudaMalloc((void**)&gpu_events_z_d, events_size);
  cudaMalloc((void**)&gpu_events_dl, events_size);

  cudaMemcpy(
      gpu_events_z_u, cpu_events_z_u, events_size, cudaMemcpyHostToDevice);
  cudaMemcpy(
      gpu_events_z_d, cpu_events_z_d, events_size, cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_events_dl, cpu_events_dl, events_size, cudaMemcpyHostToDevice);

  free(cpu_events_z_u);
  free(cpu_events_z_d);
  free(cpu_events_dl);

  for (int ib = 0; ib < n_iteration_blocks; ++ib) {
    for (int it = 0; it < n_iterations_in_block; ++it) {
      progress_callback(ib * n_iterations_in_block + it, context);

      cudaEvent_t start, stop, start_mem_time, stop_mem_time;
      float time;
      float time_all;
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      cudaEventCreate(&start_mem_time);
      cudaEventCreate(&stop_mem_time);

      cudaMemset(gpu_output, 0, output_size);
      cudaMemcpy(gpu_rho, cpu_rho, image_size, cudaMemcpyHostToDevice);

      cudaEventRecord(start);
      cudaEventRecord(start_mem_time);
#if __CUDACC__
#define reconstruction reconstruction << <blocks, threads>>>
#endif
      reconstruction(detector,
                     gpu_events_z_u,
                     gpu_events_z_d,
                     gpu_events_dl,
                     n_events,
                     gpu_output,
                     gpu_rho,
                     TEX_VAL(tex_sensitivity),
                     n_blocks,
                     n_threads_per_block);

      cudaThreadSynchronize();
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);

      cudaEventElapsedTime(&time, start, stop);
      printf("Time for the kernel: %f ms\n", time);

      // grab output
      cudaMemcpy(cpu_output, gpu_output, output_size, cudaMemcpyDeviceToHost);

      // merge image output from all blocks
      memset(cpu_rho, 0, image_size);
      for (int block = 0; block < n_blocks; ++block) {
        for (int p = 0; p < detector.total_n_pixels; ++p) {
          cpu_rho[p] += cpu_output[block * detector.total_n_pixels + p];
        }
      }

      cudaEventRecord(stop_mem_time);
      cudaEventSynchronize(stop_mem_time);
      cudaEventElapsedTime(&time_all, start_mem_time, stop_mem_time);

      printf("Time for the kernel + memory operations: %f ms\n", time_all);
    }

    output_callback(detector, ib * n_iterations_in_block, cpu_rho, context);
  }

  progress_callback(n_iteration_blocks * n_iterations_in_block, context);

#if USE_TEXTURE_OBJECT
  cudaDestroyTextureObject(tex_sensitivity);
#elif USE_TEXTURE
  cudaUnbindTexture(&tex_sensitivity);
#endif
  cudaFree(gpu_events_z_u);
  cudaFree(gpu_events_z_d);
  cudaFree(gpu_events_dl);
  cudaFree(gpu_output);
  cudaFree(gpu_rho);
  cudaFree(gpu_sensitivity);
  free(cpu_output);
  free(cpu_rho);
}

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
    int n_threads_per_block);
