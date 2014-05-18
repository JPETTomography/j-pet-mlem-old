#include <cuda_runtime.h>
#include <stdio.h>

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
                                 int event_size,
                                 int iteration_chunk,
                                 float* image_output) {

  cudaSetDevice(0);

  printf("Data Size: %d \n", event_size);
  dim3 blocks(cfg.number_of_blocks);
  dim3 threads(cfg.number_of_threads_per_block);

  cuda_kernel_config(cfg.number_of_blocks, cfg.number_of_threads_per_block);

  size_t image_sz = cfg.n_pixels * cfg.n_pixels * sizeof(float);

  float* cpu_image_buffor = (float*)malloc(image_sz * cfg.number_of_blocks);

  float* cpu_image_rho = (float*)malloc(image_sz);

  float cpu_image_sensitivity[image_sz];

  for (int i = 0; i < cfg.n_pixels * cfg.n_pixels; ++i) {

    cpu_image_rho[i] = 100.0f;
  }

  for (int i = 0; i < cfg.number_of_blocks * cfg.n_pixels * cfg.n_pixels; ++i) {

    cpu_image_buffor[i] = 0.f;
  }

  for (int px = 0; px < cfg.n_pixels; ++px) {
    for (int py = 0; py < cfg.n_pixels; ++py) {

      float2 pixel_coordiantes = pixel_center(px,
                                              py,
                                              cfg.pixel_size,
                                              cfg.pixel_size,
                                              cfg.grid_size_y_,
                                              cfg.grid_size_z_);

      cpu_image_sensitivity[px * cfg.n_pixels + py] =
          sensitivity(pixel_coordiantes.x,
                      pixel_coordiantes.y,
                      cfg.R_distance,
                      cfg.Scentilator_length / 2.0f);
    }
  }

  printf("CPU %f\n",
         sensitivity(100, 100, cfg.R_distance, cfg.Scentilator_length / 2.0f));

  float* gpu_image_buffor;
  float* gpu_image_rho;
  event<float>* gpu_event_list;
  soa_event<float>* gpu_soa_event_list;

  soa_event<float>* cpu_soa_event_list;

  cpu_soa_event_list = (soa_event<float>*)malloc(sizeof(soa_event<float>));

  cpu_soa_event_list->set_data(event_list,event_size);
  printf("%f %f %f\n",cpu_soa_event_list->z_u[0],cpu_soa_event_list->z_d[0],cpu_soa_event_list->dl[0]);
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
  cuda(Malloc, (void**)&gpu_event_list, event_size * sizeof(event<float>));

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
       event_size * sizeof(event<float>),
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

  reconstruction_2d_strip_cuda<float> << <blocks, threads>>>
      (cfg,gpu_soa_event_list, gpu_event_list, event_size, gpu_image_buffor, gpu_image_rho, tex);

  cudaThreadSynchronize();

  cudaEventRecord(stop);

  cudaEventSynchronize(stop);

  float milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);


  //33383*,16902+,1686/

  unsigned int number_of_ops_per_kernel = 51971;

  printf("Direct kernel time without memcpy %f ms\n", milliseconds);
  printf("Event list size: %d \n",event_size);
  printf("OPS per event: %d\n",number_of_ops_per_kernel);
  printf("OPS per event list: %u\n",number_of_ops_per_kernel * event_size);

  unsigned int flops_per_second = (number_of_ops_per_kernel * event_size)/(milliseconds/1000.0);

  printf("Flops_per_second: %f\n",flops_per_second/1e9);
  printf("Time in seconds: %f\n",(milliseconds/1000.0));
  printf("GFLOPS per Kernel: %f\n",(event_size * number_of_ops_per_kernel)/(milliseconds/1000.0)/1e9);

  cuda(Memcpy,
       cpu_image_buffor,
       gpu_image_buffor,
       image_sz * cfg.number_of_blocks,
       cudaMemcpyDeviceToHost);

  for (int block_id = 0; block_id < cfg.number_of_blocks; ++block_id) {
    for (int index = 0; index < cfg.n_pixels * cfg.n_pixels; ++index) {

      image_output[index] +=
          cpu_image_buffor[block_id * cfg.n_pixels * cfg.n_pixels + index];

      if (image_output[index] > 0) {
        // printf("%d %d %f\n",1,1, image_output[index]);
      }
    }
  }

  cuda(DestroyTextureObject, tex);
  cuda(Free, gpu_image_buffor);
  cuda(Free, gpu_image_rho);
  cuda(Free, texture_sensitivity_buffer);
  cuda(Free, gpu_soa_event_list);
  // free(cpu_image_sensitivity);
  free(cpu_image_buffor);
  free(cpu_image_rho);
  free(cpu_soa_event_list);
}
