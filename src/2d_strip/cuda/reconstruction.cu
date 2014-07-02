#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

#include "util/cuda/debug.h"  // catches all CUDA errors
#include "../event.h"

#include "config.h"

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
  size_t events_size = n_events * sizeof(Event<F>);

  F* cpu_output = (F*)malloc(image_size);
  F* cpu_rho = (F*)malloc(image_size);
  F* cpu_sensitivity = (F*)malloc(image_size);

  for (int i = 0; i < detector.total_n_pixels; ++i) {
    cpu_rho[i] = 100;
  }

  for (int x = 0; x < detector.n_y_pixels; ++x) {
    for (int y = 0; y < detector.n_z_pixels; ++y) {
      Point<F> point = detector.pixel_center(x, y);
      cpu_sensitivity[x * detector.n_z_pixels + y] =
          detector.sensitivity(point.x, point.y);
    }
  }

  F* gpu_output;
  F* gpu_rho;
  Event<F>* gpu_events;
  F* gpu_soa_events;
  F* cpu_soa_events;

  cpu_soa_events = (F*)malloc(events_size);
  for (int i = 0; i < n_events; ++i) {
    cpu_soa_events[i + 0 * n_events] = events[i].z_u;
    cpu_soa_events[i + 1 * n_events] = events[i].z_d;
    cpu_soa_events[i + 2 * n_events] = events[i].dl;
  }

  // declare and allocate memory
  F* gpu_sensitivity;

  size_t pitch;
  cudaMallocPitch(&gpu_sensitivity,
                  &pitch,
                  sizeof(F) * detector.n_y_pixels,
                  detector.n_z_pixels);

  cudaMemcpy2D(gpu_sensitivity,
               pitch,
               cpu_sensitivity,
               sizeof(F) * detector.n_y_pixels,
               sizeof(F) * detector.n_y_pixels,
               detector.n_z_pixels,
               cudaMemcpyHostToDevice);

  // create texture object
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
  resDesc.res.pitch2D.desc.y = 0;   // set 32 bits ONLY for float2 texture
  cudaTextureDesc texDesc;
  memset(&texDesc, 0, sizeof(texDesc));
  texDesc.readMode = cudaReadModeElementType;

  // create texture object: we only have to do this once!
  cudaTextureObject_t tex_sensitivity;
  cudaCreateTextureObject(&tex_sensitivity, &resDesc, &texDesc, NULL);

  cudaMalloc((void**)&gpu_soa_events, events_size);
  cudaMemcpy(
      gpu_soa_events, cpu_soa_events, events_size, cudaMemcpyHostToDevice);

  cudaMalloc((void**)&gpu_events, events_size);
  cudaMemcpy(gpu_events, events, events_size, cudaMemcpyHostToDevice);

  cudaMalloc((void**)&gpu_output, output_size);
  cudaMalloc((void**)&gpu_rho, image_size);

  output_callback(detector, -1, cpu_sensitivity, context);

  for (int ib = 0; ib < n_iteration_blocks; ++ib) {
    for (int it = 0; it < n_iterations_in_block; ++it) {
      progress_callback(ib * n_iterations_in_block + it, context);

      cudaMemset(gpu_output, 0, output_size);
      cudaMemcpy(gpu_rho, cpu_rho, image_size, cudaMemcpyHostToDevice);

      reconstruction_2d_strip_cuda<F> << <blocks, threads>>>
          (detector,
           gpu_soa_events,
           n_events,
           gpu_output,
           gpu_rho,
           tex_sensitivity,
           n_blocks,
           n_threads_per_block);

      cudaThreadSynchronize();

      // grab output
      cudaMemcpy(cpu_output,
                 gpu_output,
                 image_size * n_blocks,
                 cudaMemcpyDeviceToHost);

      // merge image output from all blocks
      for (int block = 0; block < n_blocks; ++block) {
        for (int p = 0; p < detector.total_n_pixels; ++p) {
          cpu_rho[p] += cpu_output[block * detector.total_n_pixels + p];
        }
      }
    }

    output_callback(detector, ib * n_iterations_in_block, cpu_rho, context);
  }

  progress_callback(n_iteration_blocks * n_iterations_in_block, context);

  cudaDestroyTextureObject(tex_sensitivity);
  cudaFree(gpu_soa_events);
  cudaFree(gpu_events);
  cudaFree(gpu_output);
  cudaFree(gpu_rho);
  cudaFree(gpu_sensitivity);
  free(cpu_soa_events);
  free(cpu_output);
  free(cpu_rho);
  free(cpu_sensitivity);
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
