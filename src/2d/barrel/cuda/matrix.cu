#include <cuda_runtime.h>
#include <stdio.h>

#include "util/cuda/debug.h"  // catches all CUDA errors
#include "util/cuda/memory.h"
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
                              float length_scale,
                              unsigned int* gpu_prng_seed,
                              int* pixel_hits) {

  bool tof = tof_step > 0;
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  util::random::tausworthe gen(&gpu_prng_seed[4 * tid]);
  util::random::uniform_real_distribution<float> one_dis(0, 1);

  __shared__ util::cuda::copy<DetectorRing> detector_ring_copier;
  detector_ring_copier = detector_ring_ptr;
  DetectorRing& detector_ring = *detector_ring_copier;

  Model model(length_scale);

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

Matrix::Matrix(const DetectorRing& detector_ring,
               int n_threads_per_block,
               int n_blocks,
               float pixel_size,
               int n_positions,
               float tof_step,
               float length_scale,
               unsigned int* prng_seed)
    : n_threads_per_block(n_threads_per_block),
      n_blocks(n_blocks),
      pixel_size(pixel_size),
      n_positions(n_positions),
      tof_step(tof_step),
      length_scale(length_scale),
      pixel_hits_count(detector_ring.n_lors * n_positions),
      pixel_hits_size(pixel_hits_count * sizeof(int)),
      output_size(n_blocks * pixel_hits_size) {

  cudaMalloc((void**)&gpu_detector_ring, sizeof(DetectorRing));
  cudaMemcpy(gpu_detector_ring,
             &detector_ring,
             sizeof(DetectorRing),
             cudaMemcpyHostToDevice);

  output = new int[n_blocks * detector_ring.n_lors * n_positions];
  cudaMalloc((void**)&gpu_output, output_size);

  int prng_seed_size = n_blocks * n_threads_per_block * 4 * sizeof(*prng_seed);
  cudaMalloc((void**)&gpu_prng_seed, prng_seed_size);
  cudaMemcpy(gpu_prng_seed, prng_seed, prng_seed_size, cudaMemcpyHostToDevice);
}

Matrix::~Matrix() {
  cudaFree(gpu_prng_seed);
  cudaFree(gpu_output);
  delete[] output;
  cudaFree(gpu_detector_ring);
}

void Matrix::operator()(Pixel pixel, int n_emissions, int* pixel_hits) {

  cudaMemset(gpu_output, 0, output_size);

#if __CUDACC__
  dim3 blocks(n_blocks);
  dim3 threads(n_threads_per_block);
#define kernel kernel << <blocks, threads>>>
#endif
  kernel(pixel,
         gpu_detector_ring,
         n_emissions,
         pixel_size,
         n_positions,
         tof_step,
         length_scale,
         gpu_prng_seed,
         gpu_output);

  cudaThreadSynchronize();
  cudaMemcpy(output, gpu_output, output_size, cudaMemcpyDeviceToHost);

  // Reduce blocks into pixel hits:
  // 1. First block can be simply copied
  memcpy(pixel_hits, output, pixel_hits_size);
  // 2. Next blocks must be reduced
  for (int block = 1; block < n_blocks; ++block) {
    for (int i = 0; i < pixel_hits_count; ++i) {
      pixel_hits[i] += output[block * pixel_hits_count + i];
    }
  }
}

}  // GPU
}  // Barrel
}  // PET2D
