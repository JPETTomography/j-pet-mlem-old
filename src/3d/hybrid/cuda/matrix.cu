#include <cuda_runtime.h>
#include <stdio.h>

#include "util/cuda/debug.h"  // catches all CUDA errors
#include "util/cuda/memory.h"
#include "util/random.h"

#include "matrix.h"

#include "3d/geometry/distribution.h"

namespace PET3D {
namespace Hybrid {
namespace GPU {

__global__ static void kernel(const float z,
                              const Pixel pixel,
                              const Scanner* scanner_ptr,
                              int n_emissions,
                              float pixel_size,
                              float length_scale,
                              unsigned int* gpu_prng_seed,
                              int* pixel_hits) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  util::random::tausworthe gen(&gpu_prng_seed[4 * tid]);
  util::random::uniform_real_distribution<float> one_dis(0, 1);
  util::random::uniform_real_distribution<float> pi_dis(0, (float)M_PI);
  Distribution::SphericalDistribution<F> direction;

  __shared__ util::cuda::copy<Scanner> scanner_shared_storage;
  scanner_shared_storage = scanner_ptr;
  Scanner& scanner = *scanner_shared_storage;

  Model model(length_scale);
  auto fov_radius2 = scanner.barrel.fov_radius() * scanner.barrel.fov_radius();

  for (int i = 0; i < n_emissions; ++i) {
    auto rx = (pixel.x + one_dis(gen)) * pixel_size;
    auto ry = (pixel.y + one_dis(gen)) * pixel_size;
    auto rz = z + one_dis(gen) * pixel_size;

    // ensure we are within a triangle
    if (rx > ry)
      continue;

    // ensure we are within FOV
    if (rx * rx + ry * ry > fov_radius2)
      continue;

    Event event(PET3D::Point<float>(rx, ry, rz), direction(gen));
    Scanner::Response response;
    auto hits = scanner.detect(gen, model, event, response);

    // do we have hit on both sides?
    if (hits >= 2) {
      auto pixel_index = response.lor.index();
      atomicAdd(&pixel_hits[pixel_index], 1);
    }
  }

  gen.save(&gpu_prng_seed[4 * tid]);
}

Matrix::Matrix(const float z,
               const Scanner& scanner,
               int n_threads_per_block,
               int n_blocks,
               float pixel_size,
               float length_scale,
               unsigned int* prng_seed)
    : z(z),
      n_threads_per_block(n_threads_per_block),
      n_blocks(n_blocks),
      pixel_size(pixel_size),
      length_scale(length_scale),
      pixel_hits_count(LOR::end_for_detectors(scanner.barrel.size()).index()),
      pixel_hits_size(pixel_hits_count * sizeof(int)) {

  cudaMalloc((void**)&gpu_scanner, sizeof(Scanner));
  cudaMemcpy(gpu_scanner, &scanner, sizeof(Scanner), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&gpu_pixel_hits, pixel_hits_size);

  int prng_seed_size = n_blocks * n_threads_per_block * 4 * sizeof(*prng_seed);
  cudaMalloc((void**)&gpu_prng_seed, prng_seed_size);
  cudaMemcpy(gpu_prng_seed, prng_seed, prng_seed_size, cudaMemcpyHostToDevice);
}

Matrix::~Matrix() {
  cudaFree(gpu_prng_seed);
  cudaFree(gpu_pixel_hits);
  cudaFree(gpu_scanner);
}

void Matrix::operator()(Pixel pixel, int n_emissions, int* pixel_hits) {

  cudaMemset(gpu_pixel_hits, 0, pixel_hits_size);

#if __CUDACC__
  dim3 blocks(n_blocks);
  dim3 threads(n_threads_per_block);
#define kernel kernel<<<blocks, threads>>>
#endif
  kernel(z,
         pixel,
         gpu_scanner,
         n_emissions,
         pixel_size,
         length_scale,
         gpu_prng_seed,
         gpu_pixel_hits);

  cudaThreadSynchronize();
  cudaMemcpy(
      pixel_hits, gpu_pixel_hits, pixel_hits_size, cudaMemcpyDeviceToHost);
}

}  // GPU
}  // Barrel
}  // PET2D
