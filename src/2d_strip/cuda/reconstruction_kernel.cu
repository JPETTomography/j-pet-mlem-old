#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

#include "../event.h"
#include "config.h"

#if EVENT_GRANULARITY
#include "reconstruction_event_granularity.cuh"
#elif WARP_GRANULARITY
#include "reconstruction_warp_granularity.cuh"
#else
#include "reconstruction_simple.cuh"
#endif

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

template <typename F>
void run_reconstruction_kernel(StripDetector<F>& detector,
                               Event<F>* events,
                               int n_events,
                               int iteration_chunk,
                               F* image_output,
                               int n_blocks,
                               int n_threads_per_block) {

  cudaSetDevice(0);

  dim3 blocks(n_blocks);
  dim3 threads(n_threads_per_block);

  size_t image_sz = detector.n_y_pixels * detector.n_z_pixels * sizeof(F);

  F* cpu_image = (F*)malloc(image_sz * n_blocks);
  F* cpu_rho = (F*)malloc(image_sz);
  F* cpu_temp_rho = (F*)malloc(image_sz);
  F cpu_sensitivity[image_sz];

  for (int i = 0; i < detector.n_y_pixels * detector.n_z_pixels; ++i) {
    cpu_rho[i] = 100;
  }

  for (int i = 0; i < n_blocks * detector.n_y_pixels * detector.n_z_pixels;
       ++i) {
    cpu_image[i] = 0;
  }

  for (int px = 0; px < detector.n_y_pixels; ++px) {
    for (int py = 0; py < detector.n_z_pixels; ++py) {
      Point<F> point = detector.pixel_center(px, py);
      cpu_sensitivity[px * detector.n_z_pixels + py] =
          detector.sensitivity(point.x, point.y);
    }
  }

  F* gpu_image;
  F* gpu_rho;
  Event<F>* gpu_events;
  SOA::Events<F>* gpu_soa_events;
  SOA::Events<F>* cpu_soa_events;

  cpu_soa_events = (SOA::Events<F>*)malloc(sizeof(SOA::Events<F>));

#ifdef OFFSET_WARP_TEST

  int offset = off;
  event<F> data_chunk[offset];

  for (int i = 0; i < offset; ++i) {
    data_chunk[i] = events[i];
  }

  cpu_soa_event_list->set_data_chunk(data_chunk, offset, event_size);
#else
  cpu_soa_events->set_data(events, n_events);
#endif
  // declare and allocate memory
  F* sensitivity_tex_buffer;

  size_t pitch;
  cudaMallocPitch(&sensitivity_tex_buffer,
                  &pitch,
                  sizeof(F) * detector.n_y_pixels,
                  detector.n_z_pixels);

  cudaMemcpy2D(sensitivity_tex_buffer,
               pitch,
               &cpu_sensitivity,
               sizeof(F) * detector.n_y_pixels,
               sizeof(F) * detector.n_y_pixels,
               detector.n_z_pixels,
               cudaMemcpyHostToDevice);

  // create texture object
  cudaResourceDesc resDesc;
  memset(&resDesc, 0, sizeof(resDesc));
  resDesc.resType = cudaResourceTypePitch2D;
  resDesc.res.pitch2D.devPtr = sensitivity_tex_buffer;
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
  cudaTextureObject_t sensitivity_tex;
  cudaCreateTextureObject(&sensitivity_tex, &resDesc, &texDesc, NULL);

  // other mallocs and allocations

  cuda(Malloc, (void**)&gpu_events, n_events * sizeof(Event<F>));
  cuda(Malloc, (void**)&gpu_image, image_sz * n_blocks);
  cuda(Malloc, (void**)&gpu_rho, image_sz);
  cuda(Malloc, (void**)&gpu_soa_events, sizeof(SOA::Events<F>));

  cuda(Memcpy,
       gpu_soa_events,
       cpu_soa_events,
       sizeof(SOA::Events<F>),
       cudaMemcpyHostToDevice);

  cuda(Memcpy,
       gpu_events,
       events,
       n_events * sizeof(Event<F>),
       cudaMemcpyHostToDevice);

  cuda(Memcpy,
       gpu_image,
       cpu_image,
       image_sz * n_blocks,
       cudaMemcpyHostToDevice);

  cuda(Memcpy, gpu_rho, cpu_rho, image_sz, cudaMemcpyHostToDevice);

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start);

  for (int i = 0; i < iteration_chunk; ++i) {

    reconstruction_2d_strip_cuda<F> << <blocks, threads>>>
        (detector,
         gpu_soa_events,
         n_events,
         gpu_image,
         gpu_rho,
         sensitivity_tex,
         n_blocks,
         n_threads_per_block);

    cudaThreadSynchronize();

    cudaEventRecord(stop);

    cudaEventSynchronize(stop);

    F milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);

    printf("Time: %f\n", milliseconds / 1000);

    cuda(Memcpy,
         cpu_image,
         gpu_image,
         image_sz * n_blocks,
         cudaMemcpyDeviceToHost);

    for (int block_id = 0; block_id < n_blocks; ++block_id) {
      for (int index = 0; index < detector.total_n_pixels; ++index) {
        image_output[i * detector.total_n_pixels + index] +=
            cpu_image[block_id * detector.total_n_pixels + index];
      }
    }

    for (int pixel = 0; pixel < n_blocks * detector.total_n_pixels; ++pixel) {
      cpu_image[pixel] = 0;
    }

    for (int pixel = 0; pixel < detector.total_n_pixels; ++pixel) {
      cpu_temp_rho[pixel] = image_output[(i * detector.total_n_pixels) + pixel];
    }

    cuda(Memcpy,
         gpu_image,
         cpu_image,
         image_sz * n_blocks,
         cudaMemcpyHostToDevice);

    cuda(Memcpy, gpu_rho, cpu_temp_rho, image_sz, cudaMemcpyHostToDevice);
  }

  // clean heap
  cuda(DestroyTextureObject, sensitivity_tex);
  cuda(Free, gpu_image);
  cuda(Free, gpu_rho);
  cuda(Free, sensitivity_tex_buffer);
  cuda(Free, gpu_soa_events);
  free(cpu_temp_rho);
  free(cpu_image);
  free(cpu_rho);
  free(cpu_soa_events);
}

template void run_reconstruction_kernel<float>(StripDetector<float>& detector,
                                               Event<float>* events,
                                               int n_events,
                                               int iteration_chunk,
                                               float* image_output,
                                               int n_blocks,
                                               int n_threads_per_block);
