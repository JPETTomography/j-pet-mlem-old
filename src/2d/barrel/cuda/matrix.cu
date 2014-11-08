#include <cuda_runtime.h>
#include <stdio.h>

#include "util/cuda/debug.h"  // catches all CUDA errors
#include "util/random.h"

#include "matrix.h"

namespace PET2D {
namespace Barrel {
namespace GPU {

__global__ static void kernel(const Pixel pixel,
                              const DetectorRing* detector_ring_ptr,
                              int n_emissions,
                              float pixel_size,
                              int n_positions,
                              float tof_step,
                              unsigned int* gpu_prng_seed,
                              int* pixel_hits) {

  bool tof = tof_step > 0;
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  util::random::tausworthe gen(&gpu_prng_seed[4 * tid]);
  util::random::uniform_real_distribution<float> one_dis(0, 1);

  __shared__
      std::aligned_storage<sizeof(DetectorRing), alignof(DetectorRing)>::type
          detector_ring_storage;
  DetectorRing& detector_ring =
      *reinterpret_cast<DetectorRing*>(&detector_ring_storage);

  if (threadIdx.x == 0) {
    memcpy(&detector_ring, detector_ring_ptr, sizeof(DetectorRing));
  }
  __syncthreads();

  Model model(0);

  for (int i = 0; i < n_emissions; ++i) {
    auto rx = (pixel.x + one_dis(gen)) * pixel_size;
    auto ry = (pixel.y + one_dis(gen)) * pixel_size;
    auto angle = one_dis(gen) * (float)M_PI;

    // ensure we are within a triangle
    if (rx > ry)
      continue;

    LOR lor;
    float position = 0;
    Event event(rx, ry, angle);
    auto hits = detector_ring.detect(gen, model, event, lor, position);

    int quantized_position = 0;
    if (tof)
      quantized_position =
          detector_ring.quantize_position(position, tof_step, n_positions);

    // do we have hit on both sides?
    if (hits >= 2) {
      auto pixel_index =
          blockIdx.x * (lor.index() * n_positions + quantized_position);
      atomicAdd(&pixel_hits[pixel_index], 1);
    }
  }

  gen.save(&gpu_prng_seed[4 * tid]);
}

bool run_matrix(Pixel pixel,
                const DetectorRing& detector_ring,
                int n_emissions,
                int n_threads_per_block,
                int n_blocks,
                float pixel_size,
                int n_positions,
                float tof_step,
                unsigned int* prng_seed,
                int* pixel_hits) {

  dim3 blocks(n_blocks);
  dim3 threads(n_threads_per_block);

  int output_size =
      n_blocks * detector_ring.n_lors * n_positions * sizeof(*pixel_hits);
  int* gpu_output;
  cudaMalloc((void**)&gpu_output, output_size);
  cudaMemset(gpu_output, 0, output_size);

  int prng_seed_size = n_blocks * n_threads_per_block * 4 * sizeof(*prng_seed);
  unsigned int* gpu_prng_seed;
  cudaMalloc((void**)&gpu_prng_seed, prng_seed_size);
  cudaMemcpy(gpu_prng_seed, prng_seed, prng_seed_size, cudaMemcpyHostToDevice);

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start);

#if __CUDACC__
#define kernel kernel << <blocks, threads>>>
#endif
  kernel(pixel,
         &detector_ring,
         n_emissions,
         pixel_size,
         n_positions,
         tof_step,
         gpu_prng_seed,
         gpu_output);

  cudaThreadSynchronize();

  cudaEventRecord(stop);
  cudaEventSynchronize(stop);

  float milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);

  cudaMemcpy(prng_seed, gpu_prng_seed, prng_seed_size, cudaMemcpyDeviceToHost);
  cudaFree(gpu_prng_seed);

  int* cpu_output = new int[n_blocks * detector_ring.n_lors * n_positions];
  cudaMemcpy(cpu_output, gpu_output, output_size, cudaMemcpyDeviceToHost);
  cudaFree(gpu_output);
  delete[] cpu_output;

  return 0;
}

}  // GPU
}  // Barrel
}  // PET2D
