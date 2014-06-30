#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

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

void gpu_reconstruction_strip_2d(CUDA::Config cfg,
                                 Event<float>* event_list,
                                 int event_size,
                                 int iteration_chunk,
                                 float* image_output,
                                 int off) {

  cudaSetDevice(0);

  dim3 blocks(cfg.number_of_blocks);
  dim3 threads(cfg.number_of_threads_per_block);

  // cuda_kernel_config(cfg.number_of_blocks,
  // cfg.number_of_threads_per_block);

  size_t image_sz = cfg.n_pixels * cfg.n_pixels * sizeof(float);

  float* cpu_image_buffor = (float*)malloc(image_sz * cfg.number_of_blocks);

  float* cpu_image_rho = (float*)malloc(image_sz);

  float* cpu_temp_rho = (float*)malloc(image_sz);

  float cpu_image_sensitivity[image_sz];

  for (int i = 0; i < cfg.n_pixels * cfg.n_pixels; ++i) {

    cpu_image_rho[i] = 100.0f;
  }

  for (int i = 0; i < cfg.number_of_blocks * cfg.n_pixels * cfg.n_pixels; ++i) {

    cpu_image_buffor[i] = 0.f;
  }

  float half_pixel_size = 0.5 * cfg.pixel_size;
  float half_grid_size = 0.5f * cfg.grid_size_y_;

  for (int px = 0; px < cfg.n_pixels; ++px) {
    for (int py = 0; py < cfg.n_pixels; ++py) {

      float2 pixel_coordiantes = pixel_center(px,
                                              py,
                                              cfg.pixel_size,
                                              cfg.pixel_size,
                                              cfg.grid_size_y_,
                                              cfg.grid_size_z_,
                                              half_grid_size,
                                              half_pixel_size);

      cpu_image_sensitivity[px * cfg.n_pixels + py] =
          sensitivity(pixel_coordiantes.x,
                      pixel_coordiantes.y,
                      cfg.R_distance,
                      cfg.Scentilator_length / 2.0f);
    }
  }

  float* gpu_image_buffor;
  float* gpu_image_rho;
  Event<float>* gpu_event_list;
  soa_event<float>* gpu_soa_event_list;

  soa_event<float>* cpu_soa_event_list;

  cpu_soa_event_list = (soa_event<float>*)malloc(sizeof(soa_event<float>));

#ifdef OFFSET_WARP_TEST

  int offset = off;

  event<float> data_chunk[offset];

  for (int i = 0; i < offset; ++i) {

    data_chunk[i] = event_list[i];
  }

  cpu_soa_event_list->set_data_chunk(data_chunk, offset, event_size);

#else

  cpu_soa_event_list->set_data(event_list, event_size);

#endif
  // declare and allocate memory
  float* texture_sensitivity_buffer;

  size_t pitch;
  cudaMallocPitch(&texture_sensitivity_buffer,
                  &pitch,
                  sizeof(float) * cfg.n_pixels,
                  cfg.n_pixels);

  cudaMemcpy2D(texture_sensitivity_buffer,
               pitch,
               &cpu_image_sensitivity,
               sizeof(float) * cfg.n_pixels,
               sizeof(float) * cfg.n_pixels,
               cfg.n_pixels,
               cudaMemcpyHostToDevice);

  // create texture object
  cudaResourceDesc resDesc;
  memset(&resDesc, 0, sizeof(resDesc));
  resDesc.resType = cudaResourceTypePitch2D;
  resDesc.res.pitch2D.devPtr = texture_sensitivity_buffer;
  resDesc.res.pitch2D.pitchInBytes = pitch;
  resDesc.res.pitch2D.width = cfg.n_pixels;
  resDesc.res.pitch2D.height = cfg.n_pixels;
  // resDesc.res.pitch2D.desc = cudaCreateChannelDesc<float>();
  resDesc.res.pitch2D.desc.f = cudaChannelFormatKindFloat;
  resDesc.res.pitch2D.desc.x = 32;  // 32 bits per channel for float texture
  resDesc.res.pitch2D.desc.y = 0;   // set 32 bits ONLY for float2 texture
  cudaTextureDesc texDesc;
  memset(&texDesc, 0, sizeof(texDesc));
  texDesc.readMode = cudaReadModeElementType;

  // create texture object: we only have to do this once!
  cudaTextureObject_t tex;
  cudaCreateTextureObject(&tex, &resDesc, &texDesc, NULL);

  // other mallocs and allocations

  cuda(Malloc, (void**)&gpu_event_list, event_size * sizeof(Event<float>));
  cuda(Malloc, (void**)&gpu_image_buffor, image_sz * cfg.number_of_blocks);
  cuda(Malloc, (void**)&gpu_image_rho, image_sz);
  cuda(Malloc, (void**)&gpu_soa_event_list, sizeof(soa_event<float>));

  cuda(Memcpy,
       gpu_soa_event_list,
       cpu_soa_event_list,
       sizeof(soa_event<float>),
       cudaMemcpyHostToDevice);

  cuda(Memcpy,
       gpu_event_list,
       event_list,
       event_size * sizeof(Event<float>),
       cudaMemcpyHostToDevice);

  cuda(Memcpy,
       gpu_image_buffor,
       cpu_image_buffor,
       image_sz * cfg.number_of_blocks,
       cudaMemcpyHostToDevice);

  cuda(Memcpy, gpu_image_rho, cpu_image_rho, image_sz, cudaMemcpyHostToDevice);

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start);

  for (int i = 0; i < iteration_chunk; ++i) {

    reconstruction_2d_strip_cuda<float> << <blocks, threads>>>
        (cfg,
         gpu_soa_event_list,
         gpu_event_list,
         event_size,
         gpu_image_buffor,
         gpu_image_rho,
         tex);

    cudaThreadSynchronize();

    cudaEventRecord(stop);

    cudaEventSynchronize(stop);

    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);

    printf("Time: %f\n", milliseconds / 1000);

    //    unsigned int number_of_ops_per_kernel = 44921;

    cuda(Memcpy,
         cpu_image_buffor,
         gpu_image_buffor,
         image_sz * cfg.number_of_blocks,
         cudaMemcpyDeviceToHost);

    for (int block_id = 0; block_id < cfg.number_of_blocks; ++block_id) {
      for (int index = 0; index < cfg.n_pixels * cfg.n_pixels; ++index) {

        image_output[(i * cfg.n_pixels * cfg.n_pixels) + index] +=
            cpu_image_buffor[block_id * cfg.n_pixels * cfg.n_pixels + index];
      }
    }

    for (int pixel = 0;
         pixel < cfg.number_of_blocks * cfg.n_pixels * cfg.n_pixels;
         ++pixel) {

      cpu_image_buffor[pixel] = 0.f;
    }

    for (int pixel = 0; pixel < cfg.n_pixels * cfg.n_pixels; ++pixel) {

      cpu_temp_rho[pixel] =
          image_output[(i * cfg.n_pixels * cfg.n_pixels) + pixel];
    }

    cuda(Memcpy,
         gpu_image_buffor,
         cpu_image_buffor,
         image_sz * cfg.number_of_blocks,
         cudaMemcpyHostToDevice);

    cuda(Memcpy, gpu_image_rho, cpu_temp_rho, image_sz, cudaMemcpyHostToDevice);
  }

  // clean heap
  cuda(DestroyTextureObject, tex);
  cuda(Free, gpu_image_buffor);
  cuda(Free, gpu_image_rho);
  cuda(Free, texture_sensitivity_buffer);
  cuda(Free, gpu_soa_event_list);
  free(cpu_temp_rho);
  free(cpu_image_buffor);
  free(cpu_image_rho);
  free(cpu_soa_event_list);
}
